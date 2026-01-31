import h5py

with h5py.File('data/inputs/integration/PBMC-23k.h5ad', 'r') as f:
    print('Root keys:', list(f.keys()))
    print()
    
    if 'obs' in f:
        print('obs:')
        if isinstance(f['obs'], h5py.Group):
            print('  Type: Group')
            print('  Keys:', list(f['obs'].keys())[:10])
            # Check for _index
            if '_index' in f['obs']:
                print('  _index shape:', f['obs']['_index'].shape)
                print('  _index dtype:', f['obs']['_index'].dtype)
        print()
    
    if 'var' in f:
        print('var:')
        if isinstance(f['var'], h5py.Group):
            print('  Type: Group')
            print('  Keys:', list(f['var'].keys())[:10])
            if '_index' in f['var']:
                print('  _index shape:', f['var']['_index'].shape)
                print('  _index dtype:', f['var']['_index'].dtype)
        print()
    
    if 'X' in f:
        print('X:')
        print('  Shape:', f['X'].shape)
        print('  Dtype:', f['X'].dtype)
