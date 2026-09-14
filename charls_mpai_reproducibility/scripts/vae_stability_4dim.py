# vae_stability_4dim.py — 10-seed stability analysis (same config as 02b)
import os, copy, numpy as np, pandas as pd, torch, torch.nn as nn, torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from scipy.linalg import orthogonal_procrustes

PROC = "/Users/a1234/Desktop/桌面/chalrs/processed"   # <-- EDIT
DEV = torch.device("cpu")
FEATS = ["tics_cog","self_health_rev","cesd_rev","physio_proxy",
         "hba1c_rev","crp_rev","cystatin_c_rev","total_chol_rev"]
HID, LAT, KLW = 64, 4, 0.5
ZCOLS = ["z%d" % (k+1) for k in range(LAT)]

class VAE(nn.Module):
    def __init__(self, d, h, l):
        super().__init__()
        self.fc1 = nn.Linear(d, h); self.fc_mu = nn.Linear(h, l)
        self.fc_logvar = nn.Linear(h, l); self.fc2 = nn.Linear(l, h); self.fc3 = nn.Linear(h, d)
    def encode(self, x):
        hh = torch.relu(self.fc1(x)); return self.fc_mu(hh), self.fc_logvar(hh)
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar); return mu + torch.randn_like(std) * std
    def decode(self, z):
        return self.fc3(torch.relu(self.fc2(z)))
    def forward(self, x):
        mu, logvar = self.encode(x); z = self.reparameterize(mu, logvar)
        return self.decode(z), mu, logvar

def loss_fn(rc, x, mu, lv, kl_w):
    recon = nn.functional.mse_loss(rc, x, reduction='sum')
    kld = -0.5 * torch.sum(1 + lv - mu.pow(2) - lv.exp())
    return recon + kl_w * kld

def train_one(X, seed, max_ep=300, patience=20):
    n = len(X); rng = np.random.RandomState(seed); perm = rng.permutation(n)
    ntr = int(n*0.8); nva = int(n*0.1)
    tr, va = perm[:ntr], perm[ntr:ntr+nva]
    torch.manual_seed(seed); np.random.seed(seed)
    model = VAE(X.shape[1], HID, LAT).to(DEV)
    opt = optim.Adam(model.parameters(), lr=1e-3)
    trl = DataLoader(TensorDataset(torch.tensor(X[tr])), batch_size=256, shuffle=True)
    best = np.inf; wait = 0; best_state = None
    for ep in range(max_ep):
        model.train()
        for (b,) in trl:
            b = b.to(DEV); opt.zero_grad()
            rc, mu, lv = model(b); loss = loss_fn(rc, b, mu, lv, KLW)
            loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            Xv = torch.tensor(X[va]).to(DEV); rc, mu, lv = model(Xv)
            vl = (nn.functional.mse_loss(rc, Xv, reduction='sum')
                  - 0.5*KLW*torch.sum(1 + lv - mu.pow(2) - lv.exp())) / len(va)
        if vl < best - 1e-4:
            best = vl; wait = 0; best_state = copy.deepcopy(model.state_dict())
        else:
            wait += 1
        if wait >= patience: break
    model.load_state_dict(best_state); return model

def loadings(X, z):
    return np.array([[np.corrcoef(X[:,j], z[:,k])[0,1] for k in range(z.shape[1])]
                     for j in range(X.shape[1])])

X1 = pd.read_csv(os.path.join(PROC, "charls_std_for_vae_m01.csv")).values.astype(np.float32)
seeds = list(range(2024, 2034))
Zs, Ls = {}, {}
for s in seeds:
    m = train_one(X1, s)
    with torch.no_grad():
        z = m.encode(torch.tensor(X1).to(DEV))[0].cpu().numpy()
    Zs[s] = z; Ls[s] = loadings(X1, z); print("seed", s, "done")

ref = Zs[seeds[0]]; Lref = Ls[seeds[0]]
phi_rows, L_all = [], []
for s in seeds:
    za = Zs[s] if s == seeds[0] else Zs[s] @ orthogonal_procrustes(Zs[s], ref)[0]
    La = loadings(X1, za); L_all.append(La)
    row = {"seed": s}
    for k in range(LAT):
        a, b = Lref[:, k], La[:, k]
        row["phi_z%d" % (k+1)] = round(float(np.sum(a*b)/np.sqrt(np.sum(a*a)*np.sum(b*b))), 3)
    phi_rows.append(row)

L_stack = np.stack(L_all)
pd.DataFrame(phi_rows).to_csv(os.path.join(PROC, "stability_phi_4dim.csv"), index=False)
pd.DataFrame(L_stack.mean(axis=0), index=FEATS, columns=ZCOLS).round(3).to_csv(
    os.path.join(PROC, "stability_mean_loadings_4dim.csv"))
pd.DataFrame(L_stack.std(axis=0), index=FEATS, columns=ZCOLS).round(3).to_csv(
    os.path.join(PROC, "stability_sd_loadings_4dim.csv"))
print(pd.DataFrame(phi_rows)); print("DONE")
