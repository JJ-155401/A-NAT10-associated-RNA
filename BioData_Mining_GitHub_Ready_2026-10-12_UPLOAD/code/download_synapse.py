import os
from pathlib import Path

import synapseclient

try:
    from _paths import EXTERNAL_DATA
except ImportError:
    # Allows direct execution from the repository's code directory.
    EXTERNAL_DATA = Path(
        os.environ.get(
            'AD_NAT10_DATA_DIR',
            Path(__file__).resolve().parents[1] / 'data' / 'external',
        )
    ).expanduser().resolve()

syn = synapseclient.Synapse()
syn.login(authToken=os.environ['SYNAPSE_AUTH_TOKEN'], silent=True)

DATA = EXTERNAL_DATA / 'downloads' / 'synapse'
DATA.mkdir(parents=True, exist_ok=True)

targets = [
    ('rosmap_tmt_protein', 'syn21971722'),   # TMT protein matrix
    ('rosmap_traits',      'syn32835854'),   # clinical traits
    ('rosmap_srm_r1',      'syn17057851'),   # SRM round 1
    ('rosmap_srm_r2',      'syn16779052'),   # SRM round 2
    ('msbb_proteomics',    'syn5759470'),    # MSBB proteomics folder
]

for name, sid in targets:
    try:
        ent = syn.get(sid, downloadLocation=str(DATA))
        fp = ent.path
        print(f'OK {name} -> {os.path.basename(fp)} ({os.path.getsize(fp)/1e6:.1f} MB)', flush=True)
    except Exception as e:
        print(f'FAIL {name} {sid}: {str(e)[:160]}', flush=True)
print('ALL_SYNAPSE_DONE', flush=True)
