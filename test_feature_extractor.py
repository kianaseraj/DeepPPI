# -*- coding: utf-8 -*-
"""
Unit tests for feature_extractor module.
Run using: pytest test_feature_extractor.py
Author: Kiana Seraj
"""

import pytest
import feature_extractor

SEQUENCE = "RKESTPHDQN"  # 10 amino acids

def test_AAC():
    aac = feature_extractor.AAC(SEQUENCE)
    assert len(aac) == 20
    assert abs(sum(aac) - 1.0) < 1e-6  # AAC should be normalized to sum to 1

def test_DPC():
    dpc = feature_extractor.DPC(SEQUENCE)
    assert len(dpc) == 400
    assert abs(sum(dpc) - 1.0) < 1e-6  # DPC should also be normalized

def test_hydrophobicity_descriptor():
    encoding = feature_extractor.hydrophobicity_descriptor(SEQUENCE)
    assert isinstance(encoding, str)
    assert set(encoding).issubset({"1", "2", "3"})
    assert len(encoding) == len(SEQUENCE)

def test_composition_descriptor():
    enc = feature_extractor.polarity_descriptor(SEQUENCE)
    comp = feature_extractor.composition_descriptor(enc)
    assert len(comp) == 3
    assert abs(sum(comp) - 1.0) < 1e-6

def test_transition_descriptor():
    enc = feature_extractor.charge_descriptor(SEQUENCE)
    trans = feature_extractor.transition_descriptor(enc)
    assert len(trans) == 3
    assert all(0.0 <= t <= 1.0 for t in trans)

def test_distribution_descriptor():
    enc = feature_extractor.secondary_structure_descriptor(SEQUENCE)
    dist = feature_extractor.distribution_descriptor(enc)
    assert len(dist) == 15
    assert all(0.0 <= d <= 1.0 for d in dist)

def test_APAAC():
    apaac = feature_extractor.APAAC(SEQUENCE)
    assert isinstance(apaac, list)
    assert len(apaac) == 80

def test_fetch_sequence(monkeypatch):
    # Mocked UniProt response for testing without internet
    class MockResponse:
        status_code = 200
        text = ">Mock\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"

    monkeypatch.setattr("requests.get", lambda url: MockResponse())
    sequence = feature_extractor.fetch_sequence("P12345")
    assert isinstance(sequence, str)
    assert sequence.startswith("MKTAY")

def test_all_property_encodings():
    descriptors = feature_extractor.physiochemical_properties
    for descriptor in descriptors:
        encoding = descriptor(SEQUENCE)
        assert isinstance(encoding, str)
        assert len(encoding) == len(SEQUENCE)
        assert set(encoding).issubset({"1", "2", "3"})
