#!/usr/bin/python

i1 = [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]] 
i2 = [[0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2], [1, 0, 2, 3], [1, 2, 3, 0], [1, 3, 0, 2], [2, 0, 1, 3], [2, 1, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2], [3, 1, 2, 0], [3, 2, 0, 1]] 
i3 = [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 3, 2, 1], [1, 0, 2, 3], [1, 0, 3, 2], [1, 2, 0, 3], [1, 2, 3, 0], [1, 3, 0, 2], [1, 3, 2, 0], [2, 0, 1, 3], [2, 0, 3, 1], [2, 1, 0, 3], [2, 1, 3, 0], [2, 3, 0, 1], [2, 3, 1, 0], [3, 0, 1, 2], [3, 0, 2, 1], [3, 1, 0, 2], [3, 1, 2, 0], [3, 2, 0, 1], [3, 2, 1, 0]] 
mnames = ['cm', 'cn', 'cr', 'cs', 'cmn', 'cmr', 'cms', 'cnm', 'cnr', 'cns', 'crm', 'crn', 'crs', 'csm', 'csn', 'csr', 'cmnr', 'cmns', 'cmrn', 'cmrs', 'cmsn', 'cmsr', 'cnmr', 'cnms', 'cnrm', 'cnrs', 'cnsm', 'cnsr', 'crmn', 'crms', 'crnm', 'crns', 'crsm', 'crsn', 'csmn', 'csmr', 'csnm', 'csnr', 'csrm', 'csrn'] 
mappings = [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 0], [1, 2], [1, 3], [2, 0], [2, 1], [2, 3], [3, 0], [3, 1], [3, 2], [0, 1, 2], [0, 1, 3], [0, 2, 1], [0, 2, 3], [0, 3, 1], [0, 3, 2], [1, 0, 2], [1, 0, 3], [1, 2, 0], [1, 2, 3], [1, 3, 0], [1, 3, 2], [2, 0, 1], [2, 0, 3], [2, 1, 0], [2, 1, 3], [2, 3, 0], [2, 3, 1], [3, 0, 1], [3, 0, 2], [3, 1, 0], [3, 1, 2], [3, 2, 0], [3, 2, 1]]

d1 = []
for seq in i1:
  plist = []
  for m in mappings:
    lm = len(m)
    pseq = [seq[pi] for pi in m]
    if lm == 1:
      for si, sm in enumerate(i1):
        if pseq == sm[:-3]:
          plist.append(si)
          break
    elif lm == 2:
      for si, sm in enumerate(i2):
        if pseq == sm[:-2]:
          plist.append(si)
          break
    else:
      for si, sm in enumerate(i3):
        if pseq == sm[:-1]:
          plist.append(si)
          break
  d1.append(plist)

d2 = []
for seq in i2:
  plist = []
  for m in mappings:
    lm = len(m)
    pseq = [seq[pi] for pi in m]
    if lm == 1:
      for si, sm in enumerate(i1):
        if pseq == sm[:-3]:
          plist.append(si)
          break
    elif lm == 2:
      for si, sm in enumerate(i2):
        if pseq == sm[:-2]:
          plist.append(si)
          break
    else:
      for si, sm in enumerate(i3):
        if pseq == sm[:-1]:
          plist.append(si)
          break
  d2.append(plist)

d3 = []
for seq in i3:
  plist = []
  for m in mappings:
    lm = len(m)
    pseq = [seq[pi] for pi in m]
    if lm == 1:
      for si, sm in enumerate(i1):
        if pseq == sm[:-3]:
          plist.append(si)
          break
    elif lm == 2:
      for si, sm in enumerate(i2):
        if pseq == sm[:-2]:
          plist.append(si)
          break
    else:
      for si, sm in enumerate(i3):
        if pseq == sm[:-1]:
          plist.append(si)
          break
  d3.append(plist)

print('typedef struct')
print('{')
for st in ['mu', 'nu', 'rho', 'sigma']:
  print('  unsigned int %s;' % st)
print('} hyp_ind_t;\n')


print('typedef struct')
print('{')
for st in mnames:
  print('  unsigned int %s;' % st)
print('} hyp_perm_t;\n')


print('static const hyp_ind_t hi1[%2d] = {' % len(i1))
for ctr, ind in enumerate(i1):
  ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
  if ctr != len(i1) - 1:
    ist += ','
  print('                                   %s' % ist) 
print('                                 };\n')


print('static const hyp_ind_t hi2[%2d] = {' % len(i2))
for ctr, ind in enumerate(i2):
  ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
  if ctr != len(i2) - 1:
    ist += ','
  print('                                   %s' % ist) 
print('                                 };\n')


print('static const hyp_ind_t hi3[%2d] = {' % len(i3))
for ctr, ind in enumerate(i3):
  ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
  if ctr != len(i3) - 1:
    ist += ','
  print('                                   %s' % ist) 
print('                                 };\n')

print('static const hyp_ind_t hp1[%2d] = {' % len(d1))
for ctr, ind in enumerate(d1):
  ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
  if ctr != len(d1) - 1:
    ist += ','
  print('                                   %s' % ist) 
print('                                 };\n')

print('static const hyp_ind_t hp2[%2d] = {' % len(d2))
for ctr, ind in enumerate(d2):
  ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
  if ctr != len(d2) - 1:
    ist += ','
  print('                                   %s' % ist) 
print('                                 };\n')

print('static const hyp_ind_t hp3[%2d] = {' % len(d3))
for ctr, ind in enumerate(d3):
  ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
  if ctr != len(d3) - 1:
    ist += ','
  print('                                   %s' % ist) 
print('                                 };\n')
