#!/bin/bash
set -e

# 1. Activate conda
source /opt/conda/etc/profile.d/conda.sh
conda activate ann

echo "=== HACGA1 Container ==="
echo "Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
echo

# If no arguments given → interactive bash
if [ $# -eq 0 ]; then
    echo "Entering interactive shell..."
    exec /bin/bash
fi

# Help
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage:"
    echo "  python /pipeline/annotation_gene/ann_gene.py [options]"
    exit 0
fi

# Execute command
exec "$@"
