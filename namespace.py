'''how I finally understood python namespaces
jan 2015'''
import namespace_import as ni
import numpy as np

print 'namespace exploration'
pi=3
print '\npi = ',pi
print 'ms.pivar = ',ni.pivar

print '\nassign to variable:'
newpi=ni.pivar
print 'newpi = ',newpi

print '\nrunning the function ms.pifun():'
ni.pifun()

print '\nms.pifun() = ',ni.pifun()

print '\nassign to variable:'
newnewpi=ni.pifun()
print 'newnewpi = ',newnewpi

print '\nnp.pi = ',np.pi
