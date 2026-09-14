# 15_dim_sensitivity.py —— 分析B：latent dimensionality 曲线 k=1,2,3,5,6
import os, copy, numpy as np, pandas as pd, torch, torch.nn as nn, torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

PROC = "/Users/a1234/Desktop/桌面/chalrs/processed"
DEV = torch.device("cpu")

class VAE(nn.Module):
    def __init__(self, d, h, l):
        super().__init__()
        self.fc1 = nn.Linear(d, h); self.fc_mu = nn.Linear(h, l)
        self.fc_logvar = nn.Linear(h, l); self.fc2 = nn.Linear(l, h); self.fc3 = nn.Linear(h, d)
    def encode(self, x):
        hh = torch.relu(self.fc1(x)); return self.fc_mu(hh), self.fc_logvar(hh)
    def reparameterize(self, mu, lv):
        return mu + torch.randn_like(torch.exp(0.5*lv)) * torch.exp(0.5*lv)
    def decode(self, z):
        return self.fc3(torch.relu(self.fc2(z)))
    def forward(self, x):
        mu, lv = self.encode(x); z = self.reparameterize(mu, lv)
        return self.decode(z), mu, lv

def train(X, lat, seed=2024, hid=64, klw=0.5, max_ep=300, patience=20):
    n = len(X); rng = np.random.RandomState(seed); perm = rng.permutation(n)
    ntr = int(n*0.8); nva = int(n*0.1)
    tr, va = perm[:ntr], perm[ntr:ntr+nva]
    torch.manual_seed(seed); np.random.seed(seed)
    m = VAE(X.shape[1], hid, lat).to(DEV)
    opt = optim.Adam(m.parameters(), lr=1e-3)
    trl = DataLoader(TensorDataset(torch.tensor(X[tr])), batch_size=256, shuffle=True)
    best = np.inf; wait = 0; bst = None
    for ep in range(max_ep):
        m.train()
        for (b,) in trl:
            b = b.to(DEV); opt.zero_grad()
            rc, mu, lv = m(b)
            loss = nn.functional.mse_loss(rc, b, reduction='sum') - 0.5*klw*torch.sum(1+lv-mu.pow(2)-lv.exp())
            loss.backward(); opt.step()
        m.eval()
        with torch.no_grad():
            Xv = torch.tensor(X[va]).to(DEV); rc, mu, lv = m(Xv)
            vl = (nn.functional.mse_loss(rc, Xv, reduction='sum')
                  - 0.5*klw*torch.sum(1+lv-mu.pow(2)-lv.exp())) / len(va)
        if vl < best - 1e-4: best = vl; wait = 0; bst = copy.deepcopy(m.state_dict())
        else: wait += 1
        if wait >= patience: break
    m.load_state_dict(bst); return m

X = pd.read_csv(os.path.join(PROC, "charls_std_for_vae_m01.csv")).values.astype(np.float32)
ids = pd.read_csv(os.path.join(PROC, "charls_ids_m01.csv"), dtype={"ID": str})
for k in [1, 2, 3, 5, 6]:
    m = train(X, k)
    with torch.no_grad():
        z = m.encode(torch.tensor(X).to(DEV))[0].cpu().numpy()
    out = pd.DataFrame(z, columns=["z%d" % (i+1) for i in range(k)])
    out["ID"] = ids["ID"].values
    out["MPAI"] = np.sqrt(np.sum(z**2, axis=1))
    out.to_csv(os.path.join(PROC, "vae_dim%d.csv" % k), index=False)
    print("dim", k, "done, val saved")
print("ALL DONE")
