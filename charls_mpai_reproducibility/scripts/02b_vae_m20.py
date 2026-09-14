# 02b_vae_m20.py — grid search, 20 VAE trainings, Procrustes alignment, Tucker phi
import os, copy, shutil, numpy as np, pandas as pd, torch, torch.nn as nn, torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from scipy.linalg import orthogonal_procrustes

PROC_DIR = "/Users/a1234/Desktop/桌面/chalrs/processed"   # <-- EDIT
DEV = torch.device("cpu")
FEATS = ["tics_cog","self_health_rev","cesd_rev","physio_proxy",
         "hba1c_rev","crp_rev","cystatin_c_rev","total_chol_rev"]

class VAE(nn.Module):
    def __init__(self, d, h=32, l=3):
        super().__init__()
        self.fc1 = nn.Linear(d, h); self.fc_mu = nn.Linear(h, l)
        self.fc_logvar = nn.Linear(h, l); self.fc2 = nn.Linear(l, h); self.fc3 = nn.Linear(h, d)
    def encode(self, x):
        h = torch.relu(self.fc1(x)); return self.fc_mu(h), self.fc_logvar(h)
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar); eps = torch.randn_like(std)
        return mu + eps * std
    def decode(self, z):
        return self.fc3(torch.relu(self.fc2(z)))
    def forward(self, x):
        mu, logvar = self.encode(x); z = self.reparameterize(mu, logvar)
        return self.decode(z), mu, logvar

def loss_fn(rc, x, mu, lv, kl_w):
    recon = nn.functional.mse_loss(rc, x, reduction='sum')
    kld = -0.5 * torch.sum(1 + lv - mu.pow(2) - lv.exp())
    return recon + kl_w * kld

def train_one(X, hid, lat, kl_w=0.5, seed=2024, max_ep=300, patience=20):
    n = len(X); rng = np.random.RandomState(seed); perm = rng.permutation(n)
    ntr = int(n*0.8); nva = int(n*0.1)
    tr, va = perm[:ntr], perm[ntr:ntr+nva]
    torch.manual_seed(seed); np.random.seed(seed)
    model = VAE(X.shape[1], hid, lat).to(DEV)
    opt = optim.Adam(model.parameters(), lr=1e-3)
    trl = DataLoader(TensorDataset(torch.tensor(X[tr])), batch_size=256, shuffle=True)
    best = np.inf; wait = 0; best_state = None; curve = []
    for ep in range(max_ep):
        model.train(); tot = 0.0
        for (b,) in trl:
            b = b.to(DEV); opt.zero_grad()
            rc, mu, lv = model(b); loss = loss_fn(rc, b, mu, lv, kl_w)
            loss.backward(); opt.step(); tot += loss.item()
        model.eval()
        with torch.no_grad():
            Xv = torch.tensor(X[va]).to(DEV); rc, mu, lv = model(Xv)
            vloss = (nn.functional.mse_loss(rc, Xv, reduction='sum')
                     - 0.5*kl_w*torch.sum(1 + lv - mu.pow(2) - lv.exp())) / len(va)
        curve.append((ep, tot/len(tr), float(vloss)))
        if vloss < best - 1e-4:
            best = vloss; wait = 0; best_state = copy.deepcopy(model.state_dict())
        else:
            wait += 1
        if wait >= patience: break
    model.load_state_dict(best_state)
    return model, np.array(curve), best

def load_m(i):
    X = pd.read_csv(os.path.join(PROC_DIR, "charls_std_for_vae_m%02d.csv" % i)).values.astype(np.float32)
    ids = pd.read_csv(os.path.join(PROC_DIR, "charls_ids_m%02d.csv" % i), dtype={"ID": str})
    assert ids["ID"].str.len().eq(12).all(), "ID error m%d" % i
    return X, ids

def loadings(X, z):
    return np.array([[np.corrcoef(X[:,j], z[:,k])[0,1] for k in range(z.shape[1])]
                     for j in range(X.shape[1])])

# ---- Step 1: grid search on imputation 01 ----
X1, _ = load_m(1)
grid = []
for lat in (2, 3, 4):
    for hid in (16, 32, 64):
        _, _, bv = train_one(X1, hid, lat)
        grid.append({"latent": lat, "hidden": hid, "val_elbo": round(float(bv), 2)})
        print("lat=%d hid=%d val_elbo=%.2f" % (lat, hid, bv))
g = pd.DataFrame(grid).sort_values("val_elbo")
g.to_csv(os.path.join(PROC_DIR, "vae_grid_search.csv"), index=False)
BEST_LAT = int(g.iloc[0]["latent"]); BEST_HID = int(g.iloc[0]["hidden"])
print("BEST: latent=%d hidden=%d" % (BEST_LAT, BEST_HID))
ZCOLS = ["z%d" % (k+1) for k in range(BEST_LAT)]

# ---- Step 2: train all 20 ----
scores = {}
for i in range(1, 21):
    X, ids = load_m(i)
    model, curve, _ = train_one(X, BEST_HID, BEST_LAT)
    np.savetxt(os.path.join(PROC_DIR, "vae_curve_m%02d.csv" % i), curve, delimiter=",")
    with torch.no_grad():
        z = model.encode(torch.tensor(X).to(DEV))[0].cpu().numpy()
    raw = pd.DataFrame(z, columns=ZCOLS)
    raw["ID"] = ids["ID"].values
    raw["MPAI"] = np.sqrt(np.sum(z**2, axis=1))
    raw.to_csv(os.path.join(PROC_DIR, "vae_scores_m%02d.csv" % i), index=False)
    scores[i] = z
    print("trained m%02d" % i)

L1 = loadings(X1, scores[1])
ld = pd.DataFrame(L1, index=FEATS, columns=ZCOLS)
print("\n=== m01 loadings ==="); print(ld.round(3))
ld.to_csv(os.path.join(PROC_DIR, "loadings_m01.csv"))

# ---- Step 3: Procrustes alignment + Tucker phi ----
phi_rows = []
for i in range(1, 21):
    z = scores[i]
    za = z if i == 1 else z @ orthogonal_procrustes(z, scores[1])[0]
    out = pd.DataFrame(za, columns=ZCOLS)
    ids = pd.read_csv(os.path.join(PROC_DIR, "charls_ids_m%02d.csv" % i), dtype={"ID": str})
    out["ID"] = ids["ID"].values
    out["MPAI"] = np.sqrt(np.sum(za**2, axis=1))
    out.to_csv(os.path.join(PROC_DIR, "vae_scores_aligned_m%02d.csv" % i), index=False)
    La = loadings(X1 if i == 1 else load_m(i)[0], za)
    row = {"m": i}
    for k in range(BEST_LAT):
        a, b = L1[:, k], La[:, k]
        row["phi_z%d" % (k+1)] = round(float(np.sum(a*b)/np.sqrt(np.sum(a*a)*np.sum(b*b))), 3)
    phi_rows.append(row)

pd.DataFrame(phi_rows).to_csv(os.path.join(PROC_DIR, "tucker_phi.csv"), index=False)
print(pd.DataFrame(phi_rows))
shutil.copy(os.path.join(PROC_DIR, "vae_scores_aligned_m01.csv"),
            os.path.join(PROC_DIR, "vae_scores.csv"))
print("DONE")
