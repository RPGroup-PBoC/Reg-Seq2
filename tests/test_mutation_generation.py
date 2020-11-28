import pytest
import wgregseq

def test_rand_seq():
    # Test input types
    with pytest.raises(TypeError):
        wgregseq.gen_rand_seq(1.5)
    with pytest.raises(TypeError):
        wgregseq.gen_rand_seq("A")

    # Test sequence length   
    seq = wgregseq.gen_rand_seq(10)
    assert len(seq) == 10

    
    