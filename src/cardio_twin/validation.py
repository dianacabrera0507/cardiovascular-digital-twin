"""Validation helpers that preserve patient-level separation."""
from __future__ import annotations
from sklearn.model_selection import GroupShuffleSplit


def patient_level_split(patient_ids, test_fraction=0.2, random_seed=2026):
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_fraction, random_state=random_seed)
    indices = list(range(len(patient_ids)))
    train, test = next(splitter.split(indices, groups=patient_ids))
    return train, test
