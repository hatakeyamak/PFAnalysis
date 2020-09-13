import matplotlib.pyplot as plt
import uproot
import numpy

def stepx(xs):
    return numpy.tile(xs, (2,1)).T.flatten()[1:-1]
def stepy(ys):
    return numpy.tile(ys, (2,1)).T.flatten()

file  = uproot.open("myfile_test.root")
file2 = uproot.open("myfile_ref.root")
file.keys()
file.allclasses()
file['ElectronEta'].show()
h = file["ElectronEta"]
h2 = file["ElectronEta_Fake"]
h_ref = file2["ElectronEta"]
h2_ref = file2["ElectronEta_Fake"]

bin_edges = h.edges # Use h.alledges if you want to include underflow and overflow bins
bin_counts = h.values # Use h.allvalues if you want to include underflow and overflow bins

bin_counts_2 = h2.values # Use h.allvalues if you want to include underflow and overflow bins
bin_counts_ref = h_ref.values # Use h.allvalues if you want to include underflow and overflow bins
bin_counts_2_ref = h2_ref.values # Use h.allvalues if you want to include underflow and overflow bins

bin_ratios = h.values/h_ref.values # Use h.allvalues if you want to include underflow and overflow bins
bin_ratios_2 = h2.values/h2_ref.values # Use h.allvalues if you want to include underflow and overflow bins

###########
# Simple one panel plotting
#
# plt.title('ElectronEta') 
# plt.xlabel('eta') 
# plt.ylabel('number of entries') 

# plt.plot(stepx(bin_edges), stepy(bin_counts), label='Gen matched')
# plt.plot(stepx(bin_edges), stepy(bin_counts_2), label='Gen no-match')
# plt.plot(stepx(bin_edges), stepy(bin_counts_ref), label='Gen matched (ref)')
# plt.plot(stepx(bin_edges), stepy(bin_counts_2_ref), label='Gen no-match (ref)')
# plt.legend()
# plt.show()

#If you want to draw a filled histogram, I usually use plt.fill_between:
#plt.fill_between(stepx(bin_edges), stepy(bin_counts))

#And for a stacked histogram, I also use fill_between:
#plt.fill_between(stepx(bin_edges), stepy(bin_counts_1))
#plt.fill_between(stepx(bin_edges), stepy(bin_counts_1), stepy(bin_counts_1 + bin_counts2))
###########

tPlot, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex='col')
tPlot.suptitle('SlimmedElectronEta') 

axes[0].plot(stepx(bin_edges), stepy(bin_counts), label='Gen matched')
axes[0].plot(stepx(bin_edges), stepy(bin_counts_2), label='Gen no-match')
axes[0].plot(stepx(bin_edges), stepy(bin_counts_ref), label='Gen matched (ref)')
axes[0].plot(stepx(bin_edges), stepy(bin_counts_2_ref), label='Gen no-match (ref)')
axes[0].legend()

axes[1].plot(stepx(bin_edges), stepy(bin_ratios), label='Gen matched')
axes[1].plot(stepx(bin_edges), stepy(bin_ratios_2), label='Gen no-match')
axes[1].legend()

axes[0].set(xlabel='eta', ylabel='number of entries')
axes[0].label_outer()
axes[0].set_ylim(bottom=0.)
axes[1].set_ylim(bottom=0.7)
axes[1].set(xlabel='eta', ylabel='test / ref')

plt.tight_layout()
plt.savefig('foo.png')
plt.savefig('foo.pdf')

plt.show()


