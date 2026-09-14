# 16_test_set_export.py —— 重构 m01 的 test 划分（与 02b 相同 seed/逻辑）
import numpy as np, pandas as pd
PROC = "/Users/a1234/Desktop/桌面/chalrs/processed"
X = pd.read_csv(PROC + "/charls_std_for_vae_m01.csv").values
ids = pd.read_csv(PROC + "/charls_ids_m01.csv", dtype={"ID": str})
n = len(X)
rng = np.random.RandomState(2024)
perm = rng.permutation(n)
ntr = int(n * 0.8); nva = int(n * 0.1)
te = perm[ntr + nva:]
out = pd.DataFrame({"ID": ids["ID"].values[te]})
out.to_csv(PROC + "/test_ids.csv", index=False)
print("test set n =", len(out))
