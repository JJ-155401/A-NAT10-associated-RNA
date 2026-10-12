"""Track A: Causal graph model + reframed narrative figure.
Reframes NAT10/ac4c from 'AD driver' to 'downstream sensitivity/amplifier'
within AD pathophysiology, reconciling: (1) MR null, (2) KO-direction reversal,
(3) ac4C-target mixed directions.

Causal chain (observational, GSE125583 n=289):
  AD pathology (Braak) --> [NF-kB up / MEF2C down] --> NAT10 down --> ac4C down
  --> subset of targets DOWN (ac4C activators: APP, MAPT, SORL1) &
      subset UP (ac4C repressors: BIN1, TREM2, CLU, APOE) --> synaptic vulnerability
"""
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os
HERE=os.path.dirname(os.path.abspath(__file__))
FIG=os.path.join(HERE,'..','figures'); os.makedirs(FIG, exist_ok=True)

fig, ax = plt.subplots(figsize=(11,6.5)); ax.axis('off')

def box(x,y,w,h,text,fc,ec='black',fs=9,tc='black'):
    ax.add_patch(plt.Rectangle((x,y),w,h,fc=fc,ec=ec, lw=1.2))
    ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=fs, color=tc, wrap=True)

# Layer 1: upstream AD drivers
box(0.02,0.70,0.20,0.10,'AD pathology\n(A-beta, tau/Braak)\n[genetic: APOE, etc.]','#fde0dd')
box(0.02,0.40,0.20,0.10,'Inflammation\n(NF-kB up)','#fde0dd')
box(0.02,0.13,0.20,0.10,'Transcriptional\nremodeling\n(MEF2C down)','#fde0dd')

# Layer 2: NAT10/ac4C (the hub we study)
box(0.40,0.40,0.20,0.14,'NAT10 writer\ndown in AD\n(log2FC=-0.16, p=5.7e-8)\nBraak rho=-0.26','#fff7bc','#d95f02',fs=8.5)
box(0.40,0.16,0.20,0.12,'ac4C landscape\ndown in AD\n(p=9.3e-7)','#fff7bc','#d95f02',fs=8.5)

# Layer 3: target subsets
box(0.70,0.58,0.26,0.12,'ac4C ACTIVATED targets\nDOWN in AD (APP, MAPT,\nSORL1) -> loss worsens','#c7e9c0','#1b9e77',fs=8)
box(0.70,0.30,0.26,0.12,'ac4C REPRESSED targets\nUP in AD (BIN1, TREM2,\nCLU, APOE) -> de-repressed','#fdd0a2','#d95f02',fs=8)

# Layer 4: outcome
box(0.70,0.02,0.26,0.12,'Synaptic\nvulnerability &\nAD progression','#d1d1e0')

# arrows
ax.annotate('',xy=(0.40,0.47),xytext=(0.22,0.75),arrowprops=dict(arrowstyle='->',color='gray'))
ax.annotate('',xy=(0.40,0.47),xytext=(0.22,0.45),arrowprops=dict(arrowstyle='->',color='gray'))
ax.annotate('',xy=(0.40,0.47),xytext=(0.22,0.18),arrowprops=dict(arrowstyle='->',color='gray'))
ax.annotate('',xy=(0.60,0.40),xytext=(0.40,0.26),arrowprops=dict(arrowstyle='->',color='#d95f02',lw=1.5))
ax.annotate('',xy=(0.70,0.64),xytext=(0.60,0.45),arrowprops=dict(arrowstyle='->',color='#1b9e77'))
ax.annotate('',xy=(0.70,0.40),xytext=(0.60,0.26),arrowprops=dict(arrowstyle='->',color='#d95f02'))
ax.annotate('',xy=(0.70,0.12),xytext=(0.70,0.30),arrowprops=dict(arrowstyle='->',color='gray'))
ax.annotate('',xy=(0.70,0.12),xytext=(0.85,0.42),arrowprops=dict(arrowstyle='->',color='gray'))

# annotations for evidence status
ax.text(0.40,0.02,'Causal status:', fontsize=8, weight='bold')
ax.text(0.40,0.00,'MR genetic: NULL (OR=1.02, p=0.33) -> no primary causality;\nNAT10 is a SENSITIVITY/AMPLIFIER, not initiator. KO-direction\nreversal explained by acute compensatory vs chronic AD remodeling.',
        fontsize=7, va='top', color='#555555')

ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_title('Reframed causal model: NAT10/ac4C as AD-downstream vulnerability amplifier\n(reconciling MR-null, KO-direction reversal, and mixed ac4C-target directions)',
             fontsize=10, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIG,'causal_model_reframed.png'), dpi=160, bbox_inches='tight')
plt.savefig(os.path.join(FIG,'causal_model_reframed.pdf'), bbox_inches='tight')
print('Saved causal_model_reframed.png')
print('DONE')
