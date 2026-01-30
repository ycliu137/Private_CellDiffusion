#!/usr/bin/env python3
"""
Test script to verify the pipeline structure and dependencies
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

print("\n" + "="*60)
print("Testing CellDiffusion + UniPort Pipeline Structure")
print("="*60)

# Test 1: Import celldiffusion
print("\n[Test 1] Importing celldiffusion module...")
try:
    import celldiffusion as cd
    print("  ✓ celldiffusion imported successfully")
except ImportError as e:
    print(f"  ✗ Failed to import celldiffusion: {e}")
    sys.exit(1)

# Test 2: Check celldiffusion API
print("\n[Test 2] Checking celldiffusion API...")
try:
    assert hasattr(cd, 'encode_features'), "Missing encode_features"
    print("  ✓ encode_features function available")
    
    assert hasattr(cd, 'inte'), "Missing inte module"
    print("  ✓ inte module available")
    
    assert hasattr(cd.inte, 'build_integration_loss_adj'), "Missing build_integration_loss_adj"
    print("  ✓ build_integration_loss_adj function available")
    
    assert hasattr(cd.inte, 'build_integration_graph'), "Missing build_integration_graph"
    print("  ✓ build_integration_graph function available")
    
    assert hasattr(cd.inte, 'integration_diffusion'), "Missing integration_diffusion"
    print("  ✓ integration_diffusion function available")
    
except AssertionError as e:
    print(f"  ✗ API check failed: {e}")
    sys.exit(1)

# Test 3: Import other dependencies
print("\n[Test 3] Checking other dependencies...")
try:
    import scanpy as sc
    print("  ✓ scanpy imported")
    
    import torch
    print(f"  ✓ torch imported (version {torch.__version__})")
    print(f"    - GPU available: {torch.cuda.is_available()}")
    
    import pandas as pd
    print("  ✓ pandas imported")
    
    import numpy as np
    print("  ✓ numpy imported")
    
except ImportError as e:
    print(f"  ✗ Failed to import dependency: {e}")
    sys.exit(1)

# Test 4: Check config file
print("\n[Test 4] Checking configuration...")
try:
    import yaml
    config_path = project_root / "PL_uniPort" / "config.yaml"
    assert config_path.exists(), f"Config file not found: {config_path}"
    
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    assert 'celldiffusion' in config, "Missing celldiffusion section in config"
    print("  ✓ celldiffusion config section found")
    
    assert 'uniport' in config, "Missing uniport section in config"
    print("  ✓ uniport config section found")
    
    # Check celldiffusion parameters
    cd_cfg = config['celldiffusion']
    required_params = ['device', 'max_epoch_fae', 'lr_fae', 'latent_size_fae', 
                       'hidden_layers_fae', 'k', 'n_edges_per_node', 'k_mnn',
                       'max_epoch_dif', 'lr_dif']
    for param in required_params:
        assert param in cd_cfg, f"Missing parameter: {param}"
    print(f"  ✓ All required celldiffusion parameters present ({len(required_params)} params)")
    
    # Check uniport parameters
    up_cfg = config['uniport']
    required_up_params = ['n_hvg_common', 'n_hvg_specific', 'mode']
    for param in required_up_params:
        assert param in up_cfg, f"Missing parameter: {param}"
    print(f"  ✓ All required uniport parameters present ({len(required_up_params)} params)")
    
except Exception as e:
    print(f"  ✗ Config check failed: {e}")
    sys.exit(1)

# Test 5: Check script files
print("\n[Test 5] Checking script files...")
try:
    scripts = {
        'run_celldiffusion.py': 'CellDiffusion integration',
        'integrate_uniport.py': 'UniPort integration',
        'plot_compare_cdif_uniport.py': 'Comparison plots'
    }
    
    scripts_dir = project_root / "PL_uniPort" / "scripts"
    for script, desc in scripts.items():
        path = scripts_dir / script
        assert path.exists(), f"Script not found: {path}"
        print(f"  ✓ {script} ({desc})")
        
except Exception as e:
    print(f"  ✗ Script check failed: {e}")
    sys.exit(1)

# Test 6: Check Snakefile
print("\n[Test 6] Checking Snakefile...")
try:
    snakefile = project_root / "PL_uniPort" / "Snakefile"
    assert snakefile.exists(), f"Snakefile not found: {snakefile}"
    print(f"  ✓ Snakefile exists")
    
    with open(snakefile) as f:
        content = f.read()
        assert 'rule run_celldiffusion' in content, "Missing run_celldiffusion rule"
        print(f"  ✓ run_celldiffusion rule found")
        
        assert 'rule integrate_uniport' in content, "Missing integrate_uniport rule"
        print(f"  ✓ integrate_uniport rule found")
        
        assert 'rule plot_compare_cdif_uniport' in content, "Missing plot_compare_cdif_uniport rule"
        print(f"  ✓ plot_compare_cdif_uniport rule found")
        
except Exception as e:
    print(f"  ✗ Snakefile check failed: {e}")
    sys.exit(1)

# Summary
print("\n" + "="*60)
print("All tests passed! ✓")
print("="*60)
print("\nPipeline is ready to run. Next steps:")
print("1. Ensure input data is in: data/inputs/integration/")
print("2. Run: cd PL_uniPort && snakemake --cores 4")
print("="*60 + "\n")
