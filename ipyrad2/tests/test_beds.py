

from pathlib import Path
from ipyrad2.assembler.beds import get_fragment_beds


OUT = Path("/tmp")
BAM = Path("/home/deren/Documents/ipyrad-tests/Ama-map-Oct5/SLH_AL_0090-contemp.sorted.filtered.bam")
assert BAM.exists(), 'where is it?'

get_fragment_beds("SLH_AL_0090-contemp", BAM, 0, 4, OUT)

