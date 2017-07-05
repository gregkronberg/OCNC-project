"""
docstring
"""

from neuron import h
import numpy as np
import matplotlib.pyplot as plt
import cell 
import stims
import param
h.load_file("stdrun.hoc")

w_ampa = 4*.00018 # ampa weight (microsiemens or micro-omhos)
w_nmda = 4*.00018 # nmda weight (microsiemens or micro-omhos)
sec_idx = [5,6,7]
seg_idx = []
for a in sec_idx:
	seg_idx.append([-1])


def run():
	"""
	docstring
	"""
	global cell1, tuft_act, fig
	cell1 = cell.Cell()
	tuft_act  = cell.Syn_act(cell1.syn_a_tuft_ampa,cell1.syn_a_tuft_nmda,sec_idx,seg_idx,w_ampa,w_nmda)

	t_vec = h.Vector()
	soma_vec = h.Vector()
	dend_vec = h.Vector()
	t_vec.record(h._ref_t)
	soma_vec.record(cell1.soma(0.5)._ref_v)
	dend_vec.record(cell1.dend_a_tuft[sec_idx[0]](0)._ref_v)
	h.dt = .01
	h.tstop = 100
	h.run()
	soma_arr = np.array(soma_vec)
	dend_arr = np.array(dend_vec)
	t_arr = np.array(t_vec)
	fig = plt.figure()
	plt.plot(t_arr,soma_arr,'k')
	plt.plot(t_arr,dend_arr)
	plt.show(fig)


if __name__ =="__main__":
	run()