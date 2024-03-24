from glob import glob
import numpy as np
import gguf

for infile in glob("*.npz"):
    outfile = infile[:-3] + "gguf"

    npz = np.load(infile)
    g = gguf.GGUFWriter(outfile, None)
    for name in npz:
        a = npz[name]
        if a.dtype == "<U3": continue
        g.add_tensor(name, a)

    g.write_header_to_file()
    g.write_kv_data_to_file()
    g.write_tensors_to_file()
    g.close()
