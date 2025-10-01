#!/usr/bin/env python

import numpy as np

MIN_SAMP_TRIM_COV = 5
ARR = np.array([
    list(bytes("AAACCCGGGTTT", "utf-8")),
    list(bytes("AAACCCGGGTTT", "utf-8")),
    list(bytes("AAACCCGGGTTT", "utf-8")),
    list(bytes("AAACCCGGGTTT", "utf-8")),
    list(bytes("NNNCCCGGG---", "utf-8")),
], dtype=np.uint8)


print(ARR.shape)

site_sample_covs = np.sum((ARR != 78) & (ARR != 45), axis=0)
print(site_sample_covs)

cov_sufficient = np.where(site_sample_covs >= MIN_SAMP_TRIM_COV)[0]
print(cov_sufficient)

trim_left = int(cov_sufficient[0])
trim_right = ARR.shape[1] - int(cov_sufficient[-1]) - 1
print(trim_left, trim_right)

tseqs = ARR[:, trim_left:ARR.shape[1] - trim_right]
print(tseqs)

tsite_sample_covs = site_sample_covs[trim_left:ARR.shape[1] - trim_right]
print(tsite_sample_covs)

print(np.sum(tsite_sample_covs > 3))
# print(trim_left, trim_right + 1)

# print(ARR[:, 0:ARR.shape[]])

print(np.concatenate([ARR, ARR], axis=1))