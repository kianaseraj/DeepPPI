import feature_extractor
#unittesting with pytest


SEQUENCE = "RKESTPHDQN"
def test_AAC():
    assert feature_extractor.AAC(SEQUENCE) == [0.0, 0.0, 10.0, 10.0, 0.0, 0.0, 10.0, 0.0, 10.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0]

def test_hydrophobicity_descriptor():
    assert feature_extractor.hydrophobicity_descriptor(SEQUENCE) == "1112222111"
