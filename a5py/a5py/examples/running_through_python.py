# coding: utf-8
import a5py.ascotpy.ascot5_main
M=a5py.ascotpy.ascot5_main.ascot5_main(input_filename=b'helloworld.h5',output_filename=b'helloworld_out.h5')
M.init()
M.read_input()
M.offload()

M.run_simulation()

M.gather_output()

print('marker summary:')
M.print_marker_summary()

M.write_output_h5()

M.free_c()

M.finalize()
