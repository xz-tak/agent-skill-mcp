#!/bin/bash
# Test the complete BioBridge skill workflow

echo "=========================================="
echo "BioBridge Skill Workflow Test"
echo "=========================================="
echo ""

# Step 1: Check environment
echo "Step 1: Checking environment..."
python scripts/ensure_env.py --quiet
if [ $? -eq 0 ]; then
    echo "✓ Environment is ready"
else
    echo "✗ Environment check failed"
    exit 1
fi
echo ""

# Step 2: Run a quick prediction test
echo "Step 2: Running prediction test..."
conda run -n biobridge python scripts/predict_link.py \
    --head IL11 \
    --head-type "gene/protein" \
    --tail-type disease \
    --topk 3 2>&1 | grep -A 5 "ENTITY MAPPING"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Prediction test successful"
else
    echo "✗ Prediction test failed"
    exit 1
fi
echo ""

# Step 3: Test characterization
echo "Step 3: Testing KG characterization..."
conda run -n biobridge python scripts/characterize_kg.py 2>&1 | head -20
if [ $? -eq 0 ]; then
    echo "✓ Characterization test successful"
else
    echo "✗ Characterization test failed"
    exit 1
fi
echo ""

echo "=========================================="
echo "✓ All workflow tests passed!"
echo "=========================================="
