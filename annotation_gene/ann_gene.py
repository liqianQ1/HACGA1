#!/bin/python3
import textwrap
import argparse
import re
import shutil
import os
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


def read_table(iput):
    key = None
    tab = []
    with open(iput, 'r') as f:
        txt = f.readline().strip()
        key = read_title(txt)
        for line in f:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            tab.append(read_value(key, line.strip()))
    return key, tab


def read_title(txt):
    key = None
    txt = txt.strip().replace(' ', '')
    if txt.startswith('#'):
        txt = txt[1:]
    key = txt.split('\t')
    return key

def read_value(key, txt):
    hash = {}
    txt = txt.strip()
    val = txt.split('\t')
    for i in range(len(key)):
        hash[key[i]] = val[i] if i < len(val) else ''
    return hash


def split_fasta_by_length(iput, S, oput):
    oput_files = []
    oput_tmp = oput.rsplit('/', 1) 
    if os.path.exists(oput + ".00"):
        oput_tmp = oput.rsplit('/', 1) 
        oput_files = [os.path.join(oput_tmp[0], file) for file in os.listdir(os.path.dirname(oput)) if file.startswith(oput.split("/")[-1] + ".")]
        return oput_files

    S = int(S or 1)

    if S == 1:
        oput_files.append(oput + ".00")
        os.symlink(iput, oput_files[0])

    else:
        fasta1 = read_fasta_all_seq(iput)
        
        T = sum(entry["len"] for entry in fasta1)
        N = math.ceil(T / S)
        i = 0
        j = 0
        k = 0
        fasta2 = []

        for entry in fasta1:
            i += 1
#            j += len(entry)
            j += entry["len"]
            fasta2.append(entry)
            if j >= N or i == len(fasta1):
                oput_files.append(os.path.join(oput_tmp[0] ,f"{oput}.{k:02d}"))
                write_fasta_all_seq(f"{oput_files[k]}", fasta2)
                fasta2 = []
                k += 1
                j = 0
    return oput_files

  
def read_fasta_all_seq(input_filename):  
    fasta_records = []  
    with open(input_filename, 'r') as input_file:  
        one_seq = ""
        for line in input_file:
            if line.startswith(">"):
                if one_seq:
                    record = read_fasta_one_seq(one_seq) 
                    fasta_records.append(record)
                    one_seq = line
                else:
                    one_seq += line
                
            else:
                one_seq += line   
        record = read_fasta_one_seq(one_seq) 
        fasta_records.append(record)
    
    return fasta_records  

def read_fasta_one_seq(one_seq):  
    lines = one_seq.split("\n") 
    header_line = lines[0] 
    seq = ""
    for line in lines:
        if line.startswith(">"):
            continue
        seq += line

    fasta_record = {  
        'id': header_line.split(" ")[0],  
        'ann':header_line.split(" ")[1] if len(header_line.split(" ")) == 2  else '',  
        'seq': seq.rstrip('*'),  
        'len': len(seq.rstrip('*'))  
    }  
       
    return fasta_record  


def write_fasta_all_seq(oput, fasta):
    with open(oput, 'w') as oput_file:
        if isinstance(fasta, list):
            for entry in fasta:
                write_fasta_one_seq(oput_file, entry)
        else:
            write_fasta_one_seq(oput_file, fasta)

def write_fasta_one_seq(oput, fasta):
    cn = 50
    if fasta['ann']:
        oput.write(f"{fasta['id']} {fasta['ann']}\n")
    else:
        oput.write(f"{fasta['id']}\n")
    seq = fasta['seq']
    for i in range(0, fasta['len'], cn):
        oput.write(seq[i:i+cn] + "\n")


def absolute_path(path, pwd=None):
    abs_path = None
    name = None
    parent_path = None

    if not pwd:
        pwd = os.getcwd()

    path = path.strip()
    pwd = pwd.strip()

    if not pwd.startswith('/'):
        raise ValueError(f"Error: {pwd} is not an absolute path")

    if not path.startswith('/'):
        path = os.path.join(pwd, path)

    old_path = path.split('/')
    new_path = []

    for part in old_path:
        if part == '' or part == '.':
            continue
        elif part == '..':
            new_path.pop()
        else:
            new_path.append(part)

    name = new_path.pop() if new_path else None
    abs_path = '/' + '/'.join(new_path) + "/*/" + name[:3] + "/" + name 
    parent_path = '/' + '/'.join(new_path) if new_path else None

    return abs_path, name, parent_path

def run_cmd(cmd):
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return (cmd, True, result.stdout)
    except subprocess.CalledProcessError as e:
        return (cmd, False, e.stderr)

        

def RNADenovo(RNAdata, fa, tran_gff):

    if os.path.exists(f"{oput_dir_protein}/01RNAdenovo/03_stringtie/stringtie.done"):
        return

    #  INDHIS_CMD
    cmd=f'''mkdir -p {oput_dir_protein}/01RNAdenovo/01_index  &&
            cd {oput_dir_protein}/01RNAdenovo/01_index/  &&
            ln -sf {fa} ref.fa  
            hisat2-build -p 8 ref.fa ./ref  &&
            mkdir -p {oput_dir_protein}/01RNAdenovo/02_hisat2/
            '''
    os.system(cmd) 

    if not os.path.exists(f"{oput_dir_protein}/01RNAdenovo/02_hisat2/hisat2.prepare.done"):
        commands = []
        with open(RNAdata)as f:
            for line in f:
                  fqsample, fq1, fq2 = line.strip().split("\t")
                  commands.append(f"hisat2 -x {oput_dir_protein}/01RNAdenovo/01_index/ref -p 4 -1 {fq1} -2 {fq2} --pen-noncansplice 1000000 |samtools view -@ 4 -bS > {oput_dir_protein}/01RNAdenovo/02_hisat2/{fqsample}.bam" )

            parallel  = int(max_parallel) // 4
     
            for i in range(0, len(commands), parallel):
                batch = commands[i:i+parallel]
                with ThreadPoolExecutor(max_workers=parallel) as executor:
                    futures = {executor.submit(run_cmd, c): c for c in batch}
                    for future in as_completed(futures):
                        cmd, success, output = future.result()
            if success:
                with open(f"{oput_dir_protein}/01RNAdenovo/02_hisat2/hisat2.prepare.done", "w"):
                    pass
            else:
                raise RuntimeError("Error executing hisat2")

    #  RNADENOVO_CMD    
    cmd=f'''mkdir -p {oput_dir_protein}/01RNAdenovo/03_stringtie  &&
            cd {oput_dir_protein}/01RNAdenovo/03_stringtie  &&
            samtools merge -@ {max_parallel} {oput_dir_protein}/01RNAdenovo/03_stringtie/merge.bam {oput_dir_protein}/01RNAdenovo/02_hisat2/*.bam  &&
            samtools sort -@ {max_parallel}  {oput_dir_protein}/01RNAdenovo/03_stringtie/merge.bam -o {oput_dir_protein}/01RNAdenovo/03_stringtie/merge.sort.bam  &&
            samtools index -@ {max_parallel} {oput_dir_protein}/01RNAdenovo/03_stringtie/merge.sort.bam   &&
            samtools flagstat -@ {max_parallel} {oput_dir_protein}/01RNAdenovo/03_stringtie/merge.sort.bam > {oput_dir_protein}/01RNAdenovo/03_stringtie/merge.flagstat.txt  &&
            stringtie -p 8 merge.sort.bam -o merge.stringtie.assemble.gtf  &&
            gffread -g {fa} merge.stringtie.assemble.gtf -w merge.stringtie.transcript.fasta  &&
            {transdecoder}/util/gtf_genome_to_cdna_fasta.pl merge.stringtie.assemble.gtf {fa} > transcripts.fasta  &&
            {transdecoder}/util/gtf_to_alignment_gff3.pl merge.stringtie.assemble.gtf > transcripts.gff3  &&
            {transdecoder}/TransDecoder.LongOrfs -t transcripts.fasta --output_dir ./   &&
            {transdecoder}/TransDecoder.Predict -t transcripts.fasta --output_dir ./   &&
            {transdecoder}/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3   &&
            sed -i "s/transdecoder/StringTie/g" transcripts.fasta.transdecoder.genome.gff3   &&
            grep -v "^#" transcripts.fasta.transdecoder.genome.gff3 > RNADenovo.gff  &&
            touch stringtie.done 
    '''
    os.system(cmd)   


def trinity_pasa(sample, fa, trinity_fa, homolog_protein, tran_gff):
    if os.path.exists(f"{oput_dir_protein}/02trinity_pasa/trinity.pasa.done"):
        return
    cmd = textwrap.dedent(f'''ln -sf /opt/conda/envs/ann/opt/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl {pasa}/scripts/Coding/cdna_alignment_orf_to_genome_orf.pl  &&
            export PATH=/opt/conda/envs/ann/bin:/opt/conda/envs/ann/opt/transdecoder/util:/opt/conda/envs/ann/opt/transdecoder:$PATH &&
            export PERL5LIB=/pipeline/annotation_gene/lib:$PERL5LIB &&
            mkdir -p {oput_dir_protein}/02trinity_pasa/ {pasa}/data &&
            cd {oput_dir_protein}/02trinity_pasa/  &&
            cp -f {root_path}/bin/alignAssembly.config  alignAssembly.config  &&
            perl  -p -i -e  's#pasa.db#{oput_dir_protein}/02trinity_pasa/pasa_{sample}#g'  alignAssembly.config  &&
            {pasa}/Launch_PASA_pipeline.pl  -g  {fa}  -t  {trinity_fa}  -c  alignAssembly.config  -C  -R  --ALIGNERS  blat,gmap,minimap2  --CPU 30  --stringent_alignment_overlap 30.0   &&
            {pasa}/scripts/pasa_asmbls_to_training_set.dbi  --pasa_transcripts_fasta  pasa_{sample}.sqlite.assemblies.fasta  --pasa_transcripts_gff3  pasa_{sample}.sqlite.pasa_assemblies.gff3  &&
            sed -i '/^$/d' pasa_{sample}.sqlite.assemblies.fasta.transdecoder.genome.gff3  &&
            {root_path}/bin/pasa_to_gff.pl  pasa_{sample}.sqlite.assemblies.fasta.transdecoder.genome.gff3  -o  pasa_{sample}.sqlite.assemblies.fasta.transdecoder.genome.gff  &&
            {root_path}/bin/gene_stats.pl  pasa_{sample}.sqlite.assemblies.fasta.transdecoder.genome.gff >stat.out  &&
            {pasa}/scripts/pasa_asmbls_to_training_set.extract_reference_orfs.pl  pasa_{sample}.sqlite.assemblies.fasta.transdecoder.genome.gff3  100  >  pasa.orfs.gff  &&
            {root_path}/bin/pasa_filter.pl pasa.orfs.gff  -o  pasa.train.gff  -cds  2  -p  pasa_{sample}.sqlite.assemblies.fasta.transdecoder.pep  -c  pasa_{sample}.sqlite.assemblies.fasta.transdecoder.cds  &&
            cp -f pasa_{sample}.sqlite.assemblies.fasta.transdecoder.genome.gff  pasa.1.end.gff  &&
            touch trinity.pasa.done
    ''')
    os.system(cmd) 

def augustus_train(sample, fa):
    if os.path.exists(f"{oput_dir_protein}/03augustus/train/augustus.train.done"):
        return
    
    gff2gbSmallDNA ,new_species, etraining, randomSplit, optimize_augustus, augustus_bin = "", "", "", "", "", augustus
    if os.path.exists(f"{augustus}/scripts"):
        gff2gbSmallDNA = f"{augustus}/scripts/gff2gbSmallDNA.pl"
        new_species = f"{augustus}/scripts/new_species.pl"
        etraining = f"{augustus}/bin/etraining"
        randomSplit = f"{augustus}/scripts/randomSplit.pl"
        optimize_augustus = f"{augustus}/scripts/optimize_augustus.pl"
        augustus_bin = f"{augustus}/bin/augustus"
    
    elif os.path.isfile(augustus):
        augustus_dir = os.path.dirname(augustus)
        gff2gbSmallDNA = f"{augustus_dir}/gff2gbSmallDNA.pl"
        new_species = f"{augustus_dir}/new_species.pl"
        etraining = f"{augustus_dir}/etraining"
        randomSplit = f"{augustus_dir}/randomSplit.pl"
        optimize_augustus = f"{augustus_dir}/optimize_augustus.pl"
        
    cmd=f'''export AUGUSTUS_CONFIG_PATH=/opt/conda/envs/ann/config  &&
            mkdir -p {oput_dir_protein}/03augustus/train  &&
            cd {oput_dir_protein}/03augustus/train  &&
            {gff2gbSmallDNA}  {oput_dir_protein}/02trinity_pasa/pasa.train.gff  {fa}  100  genes.raw.gb  &&
            {new_species}  --species={sample} --ignore  &&
            {etraining}  --species={sample}  --stopCodonExcludedFromCDS=false  genes.raw.gb  2>  validate.log  &&
            {root_path}/bin/train_filter.pl  validate.log  genes.raw.gb  -o  genes.gb  &&
            {randomSplit}  genes.gb  10  &&
            {new_species}  --species={sample} --ignore  &&
            etraining  --species={sample}  genes.gb.train  &&
            {augustus}   --species={sample}  genes.gb.test  |  tee test.out  &&
            {optimize_augustus}  --species={sample}  --cpus={max_parallel}  genes.gb.train  &&
            {etraining}  --species={sample}  genes.gb.train  &&
            augustus   --species={sample}  genes.gb.test  &&
            touch augustus.train.done
    '''
    os.system(cmd) 
    
def augustus_rna(sample, fa, fa_pre, gene_gff):
    if os.path.exists(f"{oput_dir_protein}/03augustus/augustus.done"):
        return
    
    cmd=f'''bam2hints  -intronsonly  --in={oput_dir_protein}/01RNAdenovo/03_stringtie/merge.sort.bam  --out={oput_dir_protein}/03augustus/hints.gff &&
            mkdir -p {oput_dir_protein}/03augustus/tmp
    '''
    if not os.path.exists(f"{oput_dir_protein}/03augustus/hints.gff"):
        os.system(cmd)
        
    augustus_bin = augustus
    if os.path.exists(f"{augustus}/scripts"):
        augustus_bin = f"{augustus}/bin/augustus"

    commands = [f"{augustus_bin} --species={sample}  {fa[q]}  --outfile={oput_dir_protein}/03augustus/tmp/{fa_pre[q]}.gff  --gff3=on  --uniqueGeneId=true  --hintsfile={oput_dir_protein}/03augustus/hints.gff  --alternatives-from-evidence=false  --alternatives-from-sampling=false   --allow_hinted_splicesites=atac  --extrinsicCfgFile={augustus_extrinsic}" for q in range(len(fa))]
        
    for i in range(0, len(commands), max_parallel):
        batch = commands[i:i+max_parallel]
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = {executor.submit(run_cmd, c): c for c in batch}
            for future in as_completed(futures):
                cmd, success, output = future.result()
    
    cmd=f'''export PERL5LIB=/pipeline/annotation_gene/lib:$PERL5LIB &&
            cat  {oput_dir_protein}/03augustus/tmp/*gff  >  {oput_dir_protein}/03augustus/augustus.out  &&
	       perl -i -pe 's/\ttranscript\t/\tmRNA\t/'  {oput_dir_protein}/03augustus/augustus.out  &&
	       perl {root_path}/bin/augustus_to_gff.pl  {oput_dir_protein}/03augustus/augustus.out  -o {oput_dir_protein}/03augustus/augustus.gff   &&
	       perl {root_path}/bin/gene_stats.pl {oput_dir_protein}/03augustus/augustus.gff >{oput_dir_protein}/03augustus/stat.out  &&
	       touch {oput_dir_protein}/03augustus/augustus.done
        '''
    os.system(cmd)

    
def genemarket(sample, iput1, iput2, gene_gff):
    if os.path.exists(f"{oput_dir_protein}/04genemarket/genemarket.done"):
        return
    
    cmd=f'''mkdir -p {oput_dir_protein}/04genemarket   &&
            cd {oput_dir_protein}/04genemarket   &&
            bam2hints  -intronsonly  --in={oput_dir_protein}/01RNAdenovo/03_stringtie/merge.sort.bam --out=intron.gff.tmp   &&
            {root_path}/bin/filterIntronsFindStrand.pl --genome={iput1} --introns=intron.gff.tmp --score >intron.gff   &&
            {gmes}/gmes_petap.pl  --ET intron.gff  --sequence  {iput2}  --cores  {max_parallel}  --max_gap  3000   &&
            {pasa}/misc_utilities/gtf_to_gff3_format.pl  genemark.gtf  {iput2}  >  genemark.out   &&
            {root_path}/bin/genemarkes_to_gff.pl  {iput2}  genemark.out  -o  genemarket.gff   &&
            sed -i 's/GeneMark.hmm/GeneMark.hmmET/g' genemarket.gff   &&
            {root_path}/bin/gene_stats.pl genemarket.gff > stat.out   &&
            touch genemarket.done
        '''
    os.system(cmd)
	
def genemarkep(sample, iput1, orthoDB, pro_gff):
    if os.path.exists(f"{oput_dir_protein}/05genemarkep/genemarkep.done"):
        return
    cmd=f'''export PERL5LIB=/pipeline/annotation_gene/lib:$PERL5LIB &&
            mkdir -p {oput_dir_protein}/05genemarkep   &&
            cd {oput_dir_protein}/05genemarkep  &&
            {gmes}/gmes_petap.pl --EP  --dbep {orthoDB} --sequence  {iput1} --verbose --cores={max_parallel}  &&
            {pasa}/misc_utilities/gtf_to_gff3_format.pl genemark.gtf  {iput1}  >  genemark.out   &&
            {root_path}/bin/genemarkes_to_gff.pl  {iput1}  genemark.out  -o  genemarkep.gff  &&
            {root_path}/bin/gene_stats.pl genemarkep.gff > stat.out   &&
            touch genemarkep.done
        '''
    os.system(cmd)
    

def glimmerhmm_train(sample, iput,fa, fa_pre ):
    if os.path.exists(f"{oput_dir_protein}/06glimmerhmm/train/glimmerhmm.train.done"):
        return
    
    os.makedirs(f"{oput_dir_protein}/06glimmerhmm/train/", exist_ok=True)
    cmd=f'''rm -rf {oput_dir_protein}/06glimmerhmm/train/new* '''
    os.system(cmd)
    
    commands = [f'''perl {root_path}/bin/gff_to_glimmerhmm.pl {fa[q]} {oput_dir_protein}/02trinity_pasa/pasa.train.gff  -o  {oput_dir_protein}/06glimmerhmm/train/cds.{fa_pre[q]}  && {trainGlimmerHMM} {oput_dir_protein}/06glimmerhmm/train/cds.{fa_pre[q]}.fa  {oput_dir_protein}/06glimmerhmm/train/cds.{fa_pre[q]}.lst  -i  0,40,100  -d  {oput_dir_protein}/06glimmerhmm/train/new.{fa_pre[q]}''' for q in range(len(fa))]

    for i in range(0, len(commands), int(max_parallel)):
        batch = commands[i:i+max_parallel]
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = {executor.submit(run_cmd, c): c for c in batch}
            for future in as_completed(futures):
                cmd, success, output = future.result()
    with open(f"{oput_dir_protein}/06glimmerhmm/train/glimmerhmm.train.done", 'w'):
        pass
    

def glimmerhmm_run(sample, fa, fa_pre, iput, gene_gff):
    if os.path.exists(f"{oput_dir_protein}/06glimmerhmm/glimmerhmm.done"):
        return
    os.makedirs(f"{oput_dir_protein}/06glimmerhmm/tmp/", exist_ok=True)
    commands = [f"mkdir -p {oput_dir_protein}/06glimmerhmm/tmp/{fa_pre[q]} && cd {oput_dir_protein}/06glimmerhmm/tmp/{fa_pre[q]} && perl {root_path}/bin/run_glimmerhmm.pl glimmerhmm {fa[q]} {oput_dir_protein}/06glimmerhmm/train/new.{fa_pre[q]} {oput_dir_protein}/06glimmerhmm/tmp/{fa_pre[q]}/{fa_pre[q]}.gff {oput_dir_protein}/06glimmerhmm/tmp/{fa_pre[q]}/" for q in range(len(fa))]
    for i in range(0, len(commands), int(max_parallel)):
        batch = commands[i:i+max_parallel]
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = {executor.submit(run_cmd, c): c for c in batch}
            for future in as_completed(futures):
                cmd, success, output = future.result()
                
    cmd = f'''export PERL5LIB=/pipeline/annotation_gene/lib:$PERL5LIB &&
            cat  {oput_dir_protein}/06glimmerhmm/tmp/fa*/fa*.gff  >  {oput_dir_protein}/06glimmerhmm/glimmerhmm.out  &&
            {root_path}/bin/glimmerhmm_to_gff.pl  {iput}  {oput_dir_protein}/06glimmerhmm/glimmerhmm.out  -o  {oput_dir_protein}/06glimmerhmm/glimmerhmm.gff  &&
            {root_path}/bin/gene_stats.pl {oput_dir_protein}/06glimmerhmm/glimmerhmm.gff > {oput_dir_protein}/06glimmerhmm/stat.out  &&
            touch {oput_dir_protein}/06glimmerhmm/glimmerhmm.done
        '''
    os.system(cmd)
	

def evm(sample, iput, p_gff, t_gff, g_gff):

    if os.path.exists(f"{oput_dir_protein}/07evm/evm.done"):
        return
        
    pro_gff = " ".join(["../" + file for file in p_gff])
    tran_gff = " ".join(["../" + file for file in t_gff])
    gene_gff = " ".join(["../" + file for file in g_gff])
    
    cmd=rf'''mkdir -p {oput_dir_protein}/07evm/ &&
            cd {oput_dir_protein}/07evm/ &&
            cat {pro_gff} > {oput_dir_protein}/07evm/protein_alignments.gff &&
            cat {tran_gff} > {oput_dir_protein}/07evm/transcript_alignments.gff &&
            cat {gene_gff} > {oput_dir_protein}/07evm/gene_predictions.gff &&
            cp {weights} {oput_dir_protein}/07evm/ &&
            {evidencemodeler}/EvmUtils/partition_EVM_inputs.pl  --genome  {iput}  --gene_predictions  {oput_dir_protein}/07evm/gene_predictions.gff  --protein_alignments  {oput_dir_protein}/07evm/protein_alignments.gff  --transcript_alignments  {oput_dir_protein}/07evm/transcript_alignments.gff  --segmentSize  10000000  --overlapSize  10000  --partition_listing  {oput_dir_protein}/07evm/partitions_list.out --partition_dir {oput_dir_protein}/07evm/ &&
            {evidencemodeler}/EvmUtils/write_EVM_commands.pl --genome  {iput}  --gene_predictions  {oput_dir_protein}/07evm/gene_predictions.gff  --protein_alignments  {oput_dir_protein}/07evm/protein_alignments.gff  --transcript_alignments  {oput_dir_protein}/07evm/transcript_alignments.gff  --weights {oput_dir_protein}/07evm/weights.txt  --output_file_name  evm.out  --partitions  {oput_dir_protein}/07evm/partitions_list.out  >  {oput_dir_protein}/07evm/commands.sh
        '''
    os.system(cmd)

    commands = []
    with open(f"{oput_dir_protein}/07evm/commands.sh")as f:
        for line in f:
            match = re.search(r'--exec_dir\s+(\S+)', line)
            if match:
                exec_dir = match.group(1)
                line = line.replace("-G ", f"-G {exec_dir}/").replace("-g ", f"-g {exec_dir}/").replace("-e ", f"-e {exec_dir}/").replace("-p ", f"-p {exec_dir}/")
                commands.append(line)
    for i in range(0, len(commands), max_parallel):
        batch = commands[i:i+max_parallel]
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = {executor.submit(run_cmd, c): c for c in batch}
            for future in as_completed(futures):
                cmd, success, output = future.result()
    
    cmd=f'''export PERL5LIB=/pipeline/annotation_gene/lib:$PERL5LIB &&
            cd {oput_dir_protein}/07evm/ &&
            {evidencemodeler}/EvmUtils/recombine_EVM_partial_outputs.pl  --partitions   {oput_dir_protein}/07evm/partitions_list.out  --output_file_name   {oput_dir_protein}/07evm/evm.out  &&
            {evidencemodeler}/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions   {oput_dir_protein}/07evm/partitions_list.out  --output  {oput_dir_protein}/07evm/evm.out  --genome  {iput}  &&
            cat  {oput_dir_protein}/07evm/*/evm.out.gff3  >  {oput_dir_protein}/07evm/evm.out  &&
            {root_path}/bin/evm_to_gff.pl  {oput_dir_protein}/07evm/evm.out  -o  {oput_dir_protein}/07evm/evm.gff3  &&
            touch {oput_dir_protein}/07evm/evm.done 
    '''
    os.system(cmd)
	
	
def evm_pasa(sample, iput, cufflinks_fa):

    if os.path.exists(f"{oput_dir_protein}/08evm_pasa/evm_pasa.done"):
        return
    cmd=rf'''export PERL5LIB={pasa}/SAMPLE_HOOKS:/pipeline/annotation_gene/lib:$PERL5LIB:$PERL5LIB  &&
            mkdir -p {oput_dir_protein}/08evm_pasa  &&
            cd {oput_dir_protein}/08evm_pasa  &&
            cp {root_path}/bin/alignAssembly.config  {oput_dir_protein}/08evm_pasa/alignAssembly.config  &&
            cp {root_path}/bin/pasa.annotationCompare.txt  {oput_dir_protein}/08evm_pasa/annotationCompare.config  &&
            perl  -p -i -e  's#pasa.db#{oput_dir_protein}/02trinity_pasa/pasa_{sample}#g'  alignAssembly.config  &&
	       perl  -p -i -e  's#pasa.db#{oput_dir_protein}/02trinity_pasa/pasa_{sample}#g'  annotationCompare.config  &&
	       {pasa}/scripts/Load_Current_Gene_Annotations.dbi  -c alignAssembly.config -g {iput} -P {oput_dir_protein}/07evm/evm.gff3  &&
	       {pasa}/Launch_PASA_pipeline.pl  -g  {iput}  -t  {cufflinks_fa}  -c  annotationCompare.config  -A  &&
	       {root_path}/bin/evm_to_gff.pl  *gene_structures_post_PASA_updates*.gff3  -o  {sample}.protein.gff  &&
	       {evidencemodeler}/EvmUtils/gff3_file_to_proteins.pl  {sample}.protein.gff  {iput}  prot  >  {sample}.protein.fa  &&
	       {evidencemodeler}/EvmUtils/gff3_file_to_proteins.pl  {sample}.protein.gff  {iput}  CDS   >  {sample}.cds.fa  &&
	       perl  -p -i -e  's/^(>\S+).*/{1}/'  {sample}.protein.fa  {sample}.cds.fa  &&
	       {root_path}/bin/best_gene.pl  {sample}.protein.gff  -o  {sample}.protein.best.gff  &&
	       {evidencemodeler}/EvmUtils/gff3_file_to_proteins.pl  {sample}.protein.best.gff  {iput}  prot  >  {sample}.protein.best.fa  &&
	       {evidencemodeler}/EvmUtils/gff3_file_to_proteins.pl  {sample}.protein.best.gff  {iput}  CDS   >  {sample}.cds.best.fa  &&
	       perl  -p -i -e  's/^(>\S+).*/{1}/'  {sample}.protein.best.fa  {sample}.cds.best.fa  &&
	       {root_path}/bin/gene_stats.pl {sample}.protein.best.gff >protein.best.gff.stat.out  &&
	       {root_path}/bin/gene_stats.pl {sample}.protein.gff >protein_prediction.out  &&
	       mv protein_prediction.out Genes_annotation.statistics.xls  &&
	       touch  evm_pasa.done
    '''
    os.system(cmd)


def gene_stat(sample, p_gff, t_gff, g_gff ):

    if os.path.exists(f"{oput_dir_protein}/stat.done"):
        return
    pro_gff = ", ".join(["" + file for file in p_gff])
    tran_gff = ", ".join(["" + file for file in t_gff])
    gene_gff = ", ".join(["" + file for file in g_gff])
    
    cmd=f'''cd {oput_dir_protein}  &&
    	{root_path}/bin/gff_stat.pl -g {pro_gff} > {sample}_Protein.xls  &&
	
	{root_path}/bin/gff_stat.pl -g {tran_gff} > {sample}_Transcript.xls  &&
	{root_path}/bin/gff_stat.pl -g {gene_gff} > {sample}_Abinitio.xls  &&
     touch stat.done
    '''
    os.system(cmd)

        
def get_config(config):
    config_dict = {}
    with open(config,"r") as f:
        for i in f:
            l = i.strip().split("=")
            config_dict[l[0].strip()] = l[1].strip()
    return(config_dict)
##
# trigger linguist refresh
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A pipeline for genome annotation analysis which includes denovo, homolog and the other evidence (EST, RNAseq, cDNA, et al.).")
    parser.add_argument("--Outputdir", required=True, help="Output file genome annotation analysis directory.")
    parser.add_argument("--Genome", required=True, help="Table file for sample name, genome and masked genome, forced option.")
    parser.add_argument("--RNAseq", required=True, help="Table file for sample name and fastq from transcripts, forced option.")
    parser.add_argument("--Homolog", required=True, help="Table file for homolog from relative species, forced option.")
    parser.add_argument("--max_parallel", type=int,  required=True, help="CPU number.")
    parser.add_argument("--EST", required=True, help="The result from the software Trinity and/or ISO-seq used as for EST/cDNA.")
    parser.add_argument("--weights", required=False, default="/pipeline/annotation_gene/tab/weights.txt",  help="Evidence weights used in EVM merging")
    parser.add_argument("--config", required=False, default="/pipeline/annotation_gene/config/config.txt", help="software path")
    
    args = parser.parse_args()

    Outputdir = args.Outputdir
    genome_tab = args.Genome
    homologousproteins_tab = args.Homolog
    RNAdata = args.RNAseq
    EST = args.EST
    max_parallel = args.max_parallel
    weights = args.weights
    
    config = get_config(args.config)
    transdecoder = config['transdecoder']
    pasa = config['pasa']
    augustus_extrinsic = config['augustus_extrinsic']
    augustus = config['augustus']
    gmes = config['gmes']
    pasa = config['pasa']
    trainGlimmerHMM = config['trainGlimmerHMM']
    glimmerhmm = config['glimmerhmm']
    evidencemodeler = config['evidencemodeler']
    
    
    
    oput_dir_protein = Outputdir + "/protein"
    if not os.path.exists(oput_dir_protein):
        os.makedirs(oput_dir_protein)
    RNAdata = os.path.abspath(RNAdata)
    root_path = os.path.dirname(os.path.abspath(__file__))


    keyA, sample2 = read_table(genome_tab)
    keyB, ref = read_table(homologousproteins_tab)
    keyC, transcript = read_table(EST)
    homolog_protein = ref[-1]["HOMOLOG_PROTEIN_PEP"]


    splitnA, splitnB = 100, 100

    make = []
    z = 0

    for q in range(len(sample2)):
        smp = sample2[q]
        tran = transcript[q]
        if smp["sample"] != tran["sample"]:
            raise ValueError(f"ERROR: the file {genome_tab} doesn't match {EST}")
        smp_dir = os.path.join(Outputdir, smp["sample"])
        os.makedirs(os.path.join(smp_dir, "fa1"), exist_ok=True)
        os.makedirs(os.path.join(smp_dir, "fa2"), exist_ok=True)
        
        fa1 = split_fasta_by_length(smp["scaffold"], splitnA, os.path.join(smp_dir, "fa1", "fa1"))
        fa1_pre = [os.path.basename(i) for i in fa1]
        
        fa2 = split_fasta_by_length(smp["scaffold_marked"], splitnB, os.path.join(smp_dir, "fa2", "fa2"))
        fa2_pre = [os.path.basename(i) for i in fa2]
       
        pro_gff, tran_gff, gene_gff = ["05genemarkep/genemarkep.gff"], ["01RNAdenovo/03_stringtie/RNADenovo.gff", "02trinity_pasa/pasa.1.end.gff"], ["04genemarket/genemarket.gff", "03augustus/augustus.gff", "06glimmerhmm/glimmerhmm.gff"]

        RNADenovo(RNAdata, smp["scaffold"], tran_gff)
        trinity_pasa(smp["sample"], smp["scaffold"], tran["fa"], homolog_protein, tran_gff)
        augustus_train(smp["sample"], smp["scaffold"])
        augustus_rna(smp["sample"], fa1, fa1_pre, gene_gff)
        genemarket(smp["sample"], smp["scaffold"], smp["scaffold_marked"], gene_gff)
        genemarkep(smp["sample"], smp["scaffold_marked"], homolog_protein, pro_gff)
        glimmerhmm_train(smp["sample"], smp["scaffold"], fa1, fa1_pre)
        glimmerhmm_run(smp["sample"], fa2, fa1_pre, smp["scaffold_marked"], gene_gff)
        evm(smp["sample"], smp["scaffold"], pro_gff, tran_gff, gene_gff)
        evm_pasa(smp["sample"], smp["scaffold"], tran["fa"])
        gene_stat(smp["sample"], pro_gff, tran_gff, gene_gff)
