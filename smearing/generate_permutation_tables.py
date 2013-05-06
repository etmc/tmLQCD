#!/usr/bin/python
# -*- coding: utf-8 -*-

i1 = [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]] 
i2 = [[0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2], [1, 0, 2, 3], [1, 2, 3, 0], [1, 3, 0, 2], [2, 0, 1, 3], [2, 1, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2], [3, 1, 2, 0], [3, 2, 0, 1]] 
i3 = [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 3, 2, 1], [1, 0, 2, 3], [1, 0, 3, 2], [1, 2, 0, 3], [1, 2, 3, 0], [1, 3, 0, 2], [1, 3, 2, 0], [2, 0, 1, 3], [2, 0, 3, 1], [2, 1, 0, 3], [2, 1, 3, 0], [2, 3, 0, 1], [2, 3, 1, 0], [3, 0, 1, 2], [3, 0, 2, 1], [3, 1, 0, 2], [3, 1, 2, 0], [3, 2, 0, 1], [3, 2, 1, 0]] 
mnames = ['cm', 'cn', 'cr', 'cs', 'cmn', 'cmr', 'cms', 'cnm', 'cnr', 'cns', 'crm', 'crn', 'crs', 'csm', 'csn', 'csr', 'cmnr', 'cmns', 'cmrn', 'cmrs', 'cmsn', 'cmsr', 'cnmr', 'cnms', 'cnrm', 'cnrs', 'cnsm', 'cnsr', 'crmn', 'crms', 'crnm', 'crns', 'crsm', 'crsn', 'csmn', 'csmr', 'csnm', 'csnr', 'csrm', 'csrn', 'vmap'] 
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
  seqrm = [seq[i] for i in [0, 2, 1, 3]]
  for si, sm in enumerate(i2):
    if (seq == sm) or (seqrm == sm):
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
  seqrm = [seq[i] for i in [0, 2, 1, 3]]
  for si, sm in enumerate(i2):
    if (seq == sm) or (seqrm == sm):
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
  seqrm = [seq[i] for i in [0, 2, 1, 3]]
  for si, sm in enumerate(i2):
    if (seq == sm) or (seqrm == sm):
      plist.append(si)
      break
  d3.append(plist)

with open("hex.ih", "w") as header_file:
  header_file.write('typedef struct\n')
  header_file.write('{\n')
  for st in ['mu', 'nu', 'rho', 'sigma']:
    header_file.write('  unsigned int %s;\n' % st)
  header_file.write('} hyp_ind_t;\n\n')


  header_file.write('typedef struct\n')
  header_file.write('{\n')
  for st in mnames:
    header_file.write('  unsigned int %s;\n' % st)
  header_file.write('} hyp_perm_t;\n\n')


  header_file.write('static const hyp_ind_t hi1[%2d] = {\n' % len(i1))
  for ctr, ind in enumerate(i1):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(i1) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')


  header_file.write('static const hyp_ind_t hi2[%2d] = {\n' % len(i2))
  for ctr, ind in enumerate(i2):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(i2) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')


  header_file.write('static const hyp_ind_t hi3[%2d] = {\n' % len(i3))
  for ctr, ind in enumerate(i3):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(i3) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')

  header_file.write('static const hyp_perm_t hp1[%2d] = {\n' % len(d1))
  for ctr, ind in enumerate(d1):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(d1) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')

  header_file.write('static const hyp_perm_t hp2[%2d] = {\n' % len(d2))
  for ctr, ind in enumerate(d2):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(d2) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')

  header_file.write('static const hyp_perm_t hp3[%2d] = {\n' % len(d3))
  for ctr, ind in enumerate(d3):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(d3) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')

mnames = ['cm', 'cn', 'cr', 'cmn', 'cmr', 'cnm', 'cnr', 'crm', 'crn', 'cmnr', 'cmrn', 'cnmr', 'cnrm', 'crmn', 'crnm'] 
mappings = [[1], [2], [3], [1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2], [1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

i1 = [[ 1, 2, 3], [ 2, 3, 1], [ 3, 1, 2]]
i2 = [[ 1, 2, 3], [ 1, 3, 2], [ 2, 3, 1], [ 2, 1, 3], [ 3, 1, 2], [ 3, 2, 1]]

d1 = []
for seq in i1:
  plist = []
  for m in mappings:
    lm = len(m)
    pseq = [seq[pi - 1] for pi in m]
    if lm == 1:
      for si, sm in enumerate(i1):
        if pseq == sm[:-2]:
          plist.append(si)
          break
    else:
      for si, sm in enumerate(i2):
        if pseq == sm[:-1]:
          plist.append(si)
          break
  d1.append(plist)

d2 = []
for seq in i2:
  plist = []
  for m in mappings:
    lm = len(m)
    pseq = [seq[pi - 1] for pi in m]
    if lm == 1:
      for si, sm in enumerate(i1):
        if pseq == sm[:-2]:
          plist.append(si)
          break
    else:
      for si, sm in enumerate(i2):
        if pseq == sm[:-1]:
          plist.append(si)
          break
  d2.append(plist)
  
with open("hex_3d.ih", "w") as header_file:
  header_file.write('typedef struct\n')
  header_file.write('{\n')
  for st in ['mu', 'nu', 'rho']:
    header_file.write('  unsigned int %s;\n' % st)
  header_file.write('} hyp_3d_ind_t;\n\n')


  header_file.write('typedef struct\n')
  header_file.write('{\n')
  for st in mnames:
    header_file.write('  unsigned int %s;\n' % st)
  header_file.write('} hyp_3d_perm_t;\n\n')


  header_file.write('static const hyp_3d_ind_t hi1_3d[%2d] = {\n' % len(i1))
  for ctr, ind in enumerate(i1):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(i1) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')


  header_file.write('static const hyp_3d_ind_t hi2_3d[%2d] = {\n' % len(i2))
  for ctr, ind in enumerate(i2):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(i2) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')

  header_file.write('static const hyp_3d_perm_t hp1_3d[%2d] = {\n' % len(d1))
  for ctr, ind in enumerate(d1):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(d1) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')

  header_file.write('static const hyp_3d_perm_t hp2_3d[%2d] = {\n' % len(d2))
  for ctr, ind in enumerate(d2):
    ist = '{' + ','.join(' %2d' % x for x in ind) + '}'
    if ctr != len(d2) - 1:
      ist += ','
    header_file.write('                                   %s\n' % ist) 
  header_file.write('                                 };\n\n')