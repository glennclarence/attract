import os

#TODO: fix ATTRACT-EM interface
#TODO: Symmetry help in full interface is too small

generator_version = "0.3"

parameterfiledict = {
  "ATTRACT" :  "$ATTRACTDIR/../attract.par",
  "OPLSX" : "$ATTRACTDIR/../allatom/allatom.par",
}

completion_opt = {
  ("whatif", "Protein"): "--whatif",
  ("whatif", "DNA"): "--whatif",
  ("whatif", "RNA"): "--whatif",
  ("pdb2pqr_whatif", "Protein"): "--pdb2pqr",
  ("pdb2pqr_whatif", "DNA"): "--whatif",
  ("pdb2pqr_whatif", "RNA"): "--whatif",
  #TODO: amber, cns
}  

supported_moleculetype_interactions = (
  ("Protein", "Protein"),
  ("Protein", "DNA"),
  ("Protein", "RNA"),
)  

def generate(m):     
  assert m.completion_tool not in ("amber", "cns") #TODO: AMBER and CNS completion backends not yet implemented
     
  ret = "#!/bin/bash -i\n"
  cleanupfiles = []
  if (m.header is not None):
    ret += m.header + "\n\n"
  if m.annotation is not None:
    for l in m.annotation.split("\n"):
      ret += "###" + l + "\n"
  ret += "### Generated by ATTRACT shell script generator version %s\n\n" % generator_version    
  
  ret += 'set -u -e\n'    
  ret += 'trap "kill -- -$BASHPID; $ATTRACTDIR/shm-clean" ERR EXIT\n'  
  ret += "$ATTRACTDIR/shm-clean\n\n"
  results = "result.dat result.pdb result.lrmsd result.irmsd result.fnat"
  ret += "rm -rf %s >& /dev/null\n" % results

  partnercode = ""
  filenames = []
  ensemble_lists = []
  pdbnames = []
  mappings = []
  pdbnames3 = set()
  has_tmpf = False
  
  #determine moleculetype interactions (Protein-protein, protein-RNA, protein-DNA, etc.)  
  moleculetype_interactions = set()
  for pnr,p in enumerate(m.partners):
    for pnr2,p2 in enumerate(m.partners):
      if pnr2 <= pnr: continue
      moleculetype_interaction = tuple(sorted((p.moleculetype, p2.moleculetype)))
      moleculetype_interactions.add(moleculetype_interaction)
  
  for moleculetype_interaction in moleculetype_interactions:
    if moleculetype_interaction in supported_moleculetype_interactions: continue
    if tuple(reversed(moleculetype_interaction)) in supported_moleculetype_interactions: continue
    raise ValueError("Unsupported molecule type interaction: %s - %s" % moleculetype_interaction)

  #do we need all-atom representation for RMSD calculation?
  if m.calc_irmsd or m.calc_fnat:
    aa_rmsd = True
  elif m.calc_fnat and m.rmsd_atoms == "all":
    aa_rmsd = True
  else:  
    aa_rmsd = False
  
  #do we need separate all-atom PDBs for collect, rmsd and/or iattract?     
  if m.forcefield == "OPLSX": 
    separate_aa_pdb = False
  elif m.iattract or m.collect or aa_rmsd:
    separate_aa_pdb = True
  else:
    separate_aa_pdb = False
    for p in m.partners:
      if not m.generate_modes: continue
      if m.aacontact_modes or p.moleculetype in ("DNA","RNA"):
        separate_aa_pdb = True
        break
  if separate_aa_pdb:
    aa_filenames = []
    aa_ensemble_lists = []    
  
  #do we use all-atom representation at all?
  use_aa = (m.forcefield == "OPLSX" or separate_aa_pdb == True)
  
  #do we need separate all-atom mode files?
  if not separate_aa_pdb:
    separate_aa_mode = False
  else:
    separate_aa_mode = False
    for p in m.partners:
      if not p.generate_modes: continue
      if aa_rmsd and not m.deflex:
        separate_aa_mode = True
      elif m.collect and not p.demode:
        separate_aa_mode = True
      elif m.iattract:
        separate_aa_mode = True

  modes_any = any((p.generate_modes for p in m.partners))  
  ens_any = any((p.ensemble for p in m.partners))  
   
  partnercode += """
echo '**************************************************************'
echo 'Reduce partner PDBs...'
echo '**************************************************************'
"""      
  for pnr,p in enumerate(m.partners):
    molcode = ""
    if p.moleculetype == "RNA": molcode = "--rna"
    elif p.moleculetype == "DNA": molcode = "--dna"              

    pdbname = p.pdbfile.name
    pdbname2 = os.path.split(pdbname)[1]
    pdbname3 = os.path.splitext(pdbname2)[0]    
    ensemble_list = None
    if separate_aa_pdb:
      ensemble_list_aa = None
    if not p.ensemble:
      if pdbname not in pdbnames: 
        pdbname3_0 = pdbname3
        pcount = 0
        while pdbname3 in pdbnames3:
          pcount += 1
          pdbname3 = pdbname3_0 + "-" + str(pcount)
        pdbnames3.add(pdbname3)
        pdbname4 = pdbname3 + ".pdb"
        if pdbname4 != pdbname:
          partnercode += "cat %s > %s\n" % (pdbname, pdbname4)                    
        
        #all-atom reduce; we do this even if we never use it in the docking, just to add missing atoms etc.
        pdbname_aa = pdbname3 + "-aa.pdb"
        opts = []
        if molcode: opts.append(molcode)
        if p.charged_termini: opts.append("--termini")
        if not use_aa: 
          opts.append("--heavy")
        if p.has_hydrogens:
          opts.append("--autoref")
        else:
          opts.append(completion_opt[m.completion_tool, p.moleculetype])
        opts = " ".join(opts)
        
        mapping = pdbname3 + ".mapping"
        mappings.append(mapping)
        partnercode += "python $ATTRACTDIR/../allatom/aareduce.py %s %s %s > %s\n" % (pdbname4, pdbname_aa, opts, mapping)        
        
        if m.forcefield == "ATTRACT":
          pdbname_reduced = pdbname3 + "r.pdb"
          partnercode += "python $ATTRACTTOOLS/reduce.py %s %s %s > /dev/null\n" % (pdbname_aa, pdbname_reduced, molcode)
        elif m.forcefield == "OPLSX":
          pdbname_reduced = pdbname_aa
        pdbnames.append(pdbname)
      else:
        pdbname_reduced = filenames[pdbnames.index(pdbname)]
        if separate_aa_pdb:
          pdbname_aa = aa_filenames[pdbnames.index(pdbname)]      
      filenames.append(pdbname_reduced)
      if separate_aa_pdb:
        aa_filenames.append(pdbname_aa)
    else: #p.ensemble
                  
      d = "partner%d-ensemble" % (pnr+1)
      partnercode += "rm -rf %s\n" % d
      partnercode += "mkdir %s\n" % d
      listens = d + ".list"
      dd = d + "/model"
      partnercode += "$ATTRACTTOOLS/splitmodel %s %s %d > /dev/null\n" % (p.pdbfile.name, dd, p.ensemble_size)
      partnercode += "cat /dev/null > %s\n" % listens
      if separate_aa_pdb:
        ensemble_list_aa = d + "-aa.list"
	partnercode += "cat /dev/null > %s\n" % ensemble_list_aa      
      for mnr in range(p.ensemble_size):        
        mname1 = dd + "-" + str(mnr+1) + ".pdb"
        mname2 = dd + "-" + str(mnr+1) + "r.pdb"
        mname2aa = dd + "-" + str(mnr+1) + "-aa.pdb"
        
        #all-atom reduce; we do this even if we never use it in the docking, just to add missing atoms etc.
        opts = []
        if molcode: opts.append(molcode)
        if p.charged_termini: opts.append("--termini")
        if not use_aa: 
          opts.append("--heavy")
        elif mnr > 0:  
          opts.append("--reference")
          opts.append(pdbname_reference)
        if not p.has_hydrogens:
          opts.append(completion_opt[m.completion_tool, p.moleculetype])                
        opts = " ".join(opts)
        
        mapping = "/dev/null"
        if mnr == 0:
          mapping = pdbname3 + ".mapping"
          mappings.append(mapping)        
        partnercode += "python $ATTRACTDIR/../allatom/aareduce.py %s %s %s > %s\n" % (mname1, mname2aa, opts, mapping)        
        if m.forcefield == "ATTRACT":          
          partnercode += "python $ATTRACTTOOLS/reduce.py %s %s %s > /dev/null\n" % (mname2aa, mname2, molcode)          
        elif m.forcefield == "OPLSX": 
          mname2 = mname2aa
        
        partnercode += "echo %s >> %s\n" % (mname2, listens)  
        if separate_aa_pdb:
          partnercode += "echo %s >> %s\n" % (mname2aa, ensemble_list_aa)                  
        if mnr == 0: 
          pdbname = mname2 
          pdbname_reference = pdbname
          filenames.append(pdbname)
          if separate_aa_pdb:
            pdbname_aa = mname2aa
            pdbname_reference = pdbname_aa
            aa_filenames.append(pdbname_aa)
      ensemble_list = listens      
    ensemble_lists.append(ensemble_list)
    if separate_aa_pdb:
      aa_ensemble_lists.append(ensemble_list_aa)
  partnercode += "\n"
  
  if modes_any:
    partnercode += """
echo '**************************************************************'
echo 'Assemble modes file...'
echo '**************************************************************'
cat /dev/null > hm-all.dat
""" 
    if separate_aa_mode:
      partnercode += "cat /dev/null > hm-all-aa.dat\n"
    partnercode += "\n" 
    for pnr,p in enumerate(m.partners):
      if p.generate_modes:         
        if m.forcefield == "OPLSX": 
          need_aa_modes = True         
        elif aa_rmsd and not m.deflex:
          need_aa_modes = True
        elif m.collect and not p.demode:
          need_aa_modes = True
        elif m.iattract:
          need_aa_modes = True
        else:
          need_aa_modes = False
        if need_aa_modes: assert separate_aa_pdb #if not, I made some error in the code generator logic  
        modes_file_name = "partner%d-hm.dat" % (pnr+1)
        opts = []
        moltype = ""
        if p.moleculetype == "RNA": moltype = "--rna"
        if p.moleculetype == "DNA": moltype = "--dna"
        if moltype or p.aacontact_modes:
          if moltype:
            opts.append(moltype)
          if p.aacontact_modes:
            opts.append("--aacontact")
          opts.append("--aapdb")
          opts.append(aa_filenames[pnr])
        
        opts = " ".join(opts)  
        partnercode += "python $ATTRACTTOOLS/modes.py %s %s %d > %s\n" % (filenames[pnr], opts, p.nr_modes, modes_file_name)
        if need_aa_modes:
          aa_modes_file_name = "partner%d-hm-aa.dat" % (pnr+1)
          partnercode += "python $ATTRACTTOOLS/modes.py %s %s %d > %s\n" % (aa_filenames[pnr], opts, p.nr_modes, aa_modes_file_name)
      if not p.generate_modes:
        partnercode += "echo 0 >> hm-all.dat\n"        
      else:
        partnercode += "cat %s >> hm-all.dat\n" % modes_file_name                
      if separate_aa_mode:
        if not p.generate_modes or not need_aa_modes:
          partnercode += "echo 0 >> hm-all-aa.dat\n"
        else:  
          partnercode += "cat %s >> hm-all-aa.dat\n" % aa_modes_file_name        
    partnercode += "\n"  
  
  if use_aa: 
    if not separate_aa_pdb:
      aa_filenames = filenames
    if ens_any:
      if not separate_aa_pdb:
        aa_ensemble_lists = ensemble_lists
    if modes_any:
      if separate_aa_mode:
        aa_modesfile = "hm-all-aa.dat"    
      else:
        aa_modesfile = "hm-all.dat"
    rmsd_filenames = aa_filenames
  else:
    rmsd_filenames = filenames
    
  if len(m.partners) == 2:
    partnerfiles = " ".join(filenames)
  else:
    partnerfiles = "partners.pdb"
    partnercode += """
echo '**************************************************************'
echo 'Concatenate partner PDBs...'
echo '**************************************************************'
echo %d > partners.pdb
""" % len(m.partners)
    for f in filenames:
      partnercode += "grep ATOM %s >> partners.pdb\n" % f
      partnercode += "echo TER >> partners.pdb\n"
          
  ret += """
#name of the run
name=%s 
""" % m.runname
  ffpar = parameterfiledict[m.forcefield]
  params = "\"" + ffpar + " " + partnerfiles
  scoreparams = params + " --score --fix-receptor"
  gridparams = ""
  if m.fix_receptor: params += " --fix-receptor" 
  if modes_any: 
    ps = " --modes hm-all.dat"
    params += ps
    scoreparams += ps
  gridfiles = {}
  ret_shm = ""
  for g in m.grids:
    gheader = ""
    if m.np > 1: gheader = "header"
    extension = ".tgrid" if g.torque else ".grid"
    v = g.gridname.strip() + extension + gheader
    gridfiles[g.gridname.strip()] = (v, g.torque)
  grid_used = {}
  for pnr,p in enumerate(m.partners):
    if p.gridname is not None:
      v, is_torque = gridfiles[p.gridname.strip()]
      if v in grid_used: 
        v = grid_used[v]
      else:
        grid_used[v] = pnr+1
      gridoption = "--torquegrid" if is_torque else "--grid"
      gridparams += " %s %d %s" % (gridoption, pnr+1, str(v))
    if ensemble_lists[pnr] is not None:
      ps = " --ens %d %s" % (pnr+1, ensemble_lists[pnr])
      params += ps
      scoreparams += ps
  if m.ghost:
    params += " --ghost"
  if m.ghost_ligands:
    params += " --ghost-ligands"    
  if m.gravity:
    params += " --gravity %d" % m.gravity
  if m.rstk != 0.01:
    params += " --rstk %s" % str(m.rstk)
  if m.dielec == "cdie": 
    ps = " --cdie"
    params += ps  
    scoreparams += ps  
  if m.epsilon != 15: 
    ps = " --epsilon %s" % (str(m.epsilon))
    params += ps
    scoreparams += ps  
  has_restraints = False 
  has_haddock_restraints = False
  if m.restraints_file is not None:
    ps = " --rest %s" % (str(m.restraints_file.name))
    params += ps
    has_restraints = True
  for p in m.partners:
    if p.haddock_restraints is not None:
      has_haddock_restraints = True
  if has_haddock_restraints:    
    haddock_restraints_filename = "haddock-restraints-%s.rest" % m.runname
    ps = " --rest %s" % haddock_restraints_filename
    params += ps
    has_restraints = True
  if has_restraints:
    ps_score = ps + " --restweight %s"  % (str(m.restraints_score_weight))
    scoreparams += ps_score
        
  for sym in m.symmetries:        
    if sym.symmetry_axis is None: #distance-restrained symmetry      
      symcode = len(sym.partners)
      if sym.symmetry == "Dx": symcode = -4
      partners = " ".join([str(v) for v in sym.partners])
      params += " --sym %d %s" % (symcode, partners)
    else: #generative symmetry
      symcode = sym.symmetry_fold
      partner = sym.partners[0]
      s = sym.symmetry_axis
      sym_axis = "%.3f %.3f %.3f" % (s.x, s.y, s.z)
      s = sym.symmetry_origin
      sym_origin = "%.3f %.3f %.3f" % (s.x, s.y, s.z)      
      params += " --axsym %d %d %s %s" % (partner, symcode, sym_axis, sym_origin)
  paramsprep = params.replace("--fix-receptor","").replace("--ghost-ligands","").replace("--ghost","").replace("  ", " ") + " --ghost"
  params += "\""
  paramsprep += "\""
  scoreparams += "\""  
  ret += """
#docking parameters
params=%s
paramsprep=%s
scoreparams=%s
"""  % (params, paramsprep, scoreparams)
  if len(gridparams):
    ret += """
#grid parameters
gridparams="%s"
"""  % gridparams
  
  if m.np > 1:
    if m.jobsize == 0:
      parals = "--np %d --chunks %d" % (m.np, m.np)
    else:
      parals = "--np %d --jobsize %d" % (m.np, m.jobsize)
    ret += """
#parallelization parameters
parals="%s"
"""  % (parals)
        
  ret += "if [ 1 -eq 1 ]; then ### move and change to disable parts of the protocol\n"  
  ret += partnercode  
  ret += ret_shm
  
  #determine flexibility parameters for fix_receptor and deredundant during docking
  flexpar1 = ""
  flexpar_iattract = ""
  if ens_any:
    flexpar1 += " --ens"
    if m.iattract: 
      flexpar_iattract += " --ens"
    for p in m.partners:
      ensemble_size = p.ensemble_size
      flexpar1 += " %d" % ensemble_size
      if m.iattract: 
        flexpar_iattract += " %d" % ensemble_size        
    if m.deredundant_ignorens: flexpar1 += " --ignorens"
  if modes_any:
    flexpar1 += " --modes"
    for p in m.partners:
      nr_modes = p.nr_modes
      flexpar1 += " %d" % nr_modes
  
  if has_haddock_restraints:
      ret += """
echo '**************************************************************'
echo 'Generate HADDOCK restraints...'
echo '**************************************************************'
"""
      
      air_filenames = []
      for pnr, p in enumerate(m.partners):        
        if p.haddock_restraints:
          f_act = "haddock-active-%s-partner%d.reslist" % (m.runname, pnr+1)
          actstr = "\n".join([str(r) for r in p.haddock_restraints.activereslist])
          open(f_act, "w").write(actstr)
          air_filenames.append(f_act)

          f_pass = "haddock-passive-%s-partner%d.reslist" % (m.runname, pnr+1)
          passtr = "\n".join([str(r) for r in p.haddock_restraints.passivereslist])
          open(f_pass, "w").write(passtr)
          air_filenames.append(f_pass)
        
          air_filenames.append(filenames[pnr])
          air_filenames.append(mappings[pnr])
          air_filenames.append("\\\n" + " " * 27)
          
      dist = 2.0
      if m.forcefield == "ATTRACT": dist = 3.0
      #TODO: version of air.py for multi-body docking
      ret += "python $ATTRACTTOOLS/air.py %s 0.5 %s > %s\n" % (" ".join(air_filenames), dist, haddock_restraints_filename)
    
  if m.search == "syst" or m.search == "custom":
    if m.search == "syst" or m.start_structures_file is None:
      ret += """
echo '**************************************************************'
echo 'Generate starting structures...'
echo '**************************************************************'
"""
      rotfil = "$ATTRACTDIR/../rotation.dat"
      if m.rotations_file is not None:
        rotfil = m.rotations_file.name
      if rotfil != "rotation.dat":
        ret += "cat %s > rotation.dat\n" % rotfil
      if m.translations_file is not None:
        transfil = m.translations_file.name
        if transfil != "translate.dat":
          ret += "cat %s > translate.dat\n" % transfil
      else:
        ret += "$ATTRACTDIR/translate %s %s > translate.dat\n" % \
         (filenames[0], filenames[1])
      ret += "$ATTRACTDIR/systsearch > systsearch.dat\n"
      ret += "start=systsearch.dat\n"
      start = "systsearch.dat"
    else:
      ret += """
#starting structures
start=%s 
""" % m.start_structures_file.name
      start = m.start_structures_file.name
  elif m.search == "random":
    ret += """
echo '**************************************************************'
echo 'Generate starting structures...'
echo '**************************************************************'
"""    
    fixre = ""
    if m.fix_receptor: fixre = " --fix-receptor"
    ret += "python $ATTRACTTOOLS/randsearch.py %d %d%s > randsearch.dat\n" % \
     (len(m.partners), m.structures, fixre)
    ret += "start=randsearch.dat\n"    
    start = "randsearch.dat"
  else:
    raise Exception("Unknown value")
  ret += "\n"  
  
  inp = "$start"
  if any((p.ensemblize for p in m.partners if p.ensemblize not in (None,"custom"))):
    start0 = start
    ret += """
echo '**************************************************************'
echo 'ensemble search:' 
echo ' add ensemble conformations to the starting structures'
echo '**************************************************************'
"""
    for pnr, p in enumerate(m.partners):
      if p.ensemblize in (None, "custom"): continue
      if p.ensemblize == "all":
        ret += """
echo '**************************************************************'
echo 'multiply the number of structures by the ensemble for ligand %d'
echo '**************************************************************'
""" % (pnr+1)
      elif p.ensemblize == "random":
        ret += """
echo '**************************************************************'
echo 'random ensemble conformation in ligand %d for each starting structure'
echo '**************************************************************'
""" % (pnr+1)
      else: 
        raise ValueError(p.ensemblize)
      if start == start0 and start not in ("systsearch.dat", "randsearch.dat"):
        start2 = "start-ens%d.dat" % (pnr+1)
      else:
        start2 = os.path.splitext(start)[0] + "-ens%d.dat" % (pnr+1)
      ret += "python $ATTRACTTOOLS/ensemblize.py %s %d %d %s  > %s\n" % \
       (start, p.ensemble_size, pnr+1, p.ensemblize, start2)
      start = start2
      inp = start 
    
  for g in m.grids:
    gridfile = gridfiles[g.gridname.strip()][0]
    if gridfile not in grid_used: continue
    ret += """
echo '**************************************************************'
echo 'calculate %s grid'
echo '**************************************************************'
""" % g.gridname
    partner = grid_used[gridfile]-1
    f = filenames[partner]
    f0 = os.path.splitext(f)[0]
    tomp = ""
    if g.torque: tomp += "-torque"
    if g.omp: tomp += "-omp"
    tail = ""
    if m.np > 1: 
      tail += " --shm"
    direc = ">"
    for fnr in range(len(filenames)):
      if fnr == partner: continue
      ret += "awk '{print substr($0,58,2)}' %s | sort -nu %s %s.alphabet\n" % (filenames[fnr], direc, g.gridname.strip())
      direc = ">>"
    tail += " --alphabet %s.alphabet" % g.gridname.strip()  
    if m.dielec == "cdie": tail += " --cdie"
    if m.epsilon != 15: tail += " --epsilon %s" % (str(m.epsilon))    
    if g.calc_potentials == False: tail += " -calc-potentials=0"
    ret += "$ATTRACTDIR/make-grid%s %s %s %s %s %s %s\n\n" % \
     (tomp, f, ffpar, g.plateau_distance, g.neighbour_distance, gridfile, tail)
  ret += """
echo '**************************************************************'
echo 'Docking'
echo '**************************************************************'
"""
  itparams0 = ""
  if m.atomdensitygrids is not None:
    for g in m.atomdensitygrids:
      itparams0 += " --atomdensitygrid %s %s %s" % (g.voxelsize, g.dimension, g.forceconstant)
  ordinals = ["1st", "2nd", "3rd",] + ["%dth" % n for n in range(4,51)]
  iterations = []
  outp = ""
  if m.nr_iterations is None: m.nr_iterations = len(m.iterations)
  for n in range(m.nr_iterations):  
    if m.iterations is None or len(m.iterations) <= n:
      it = AttractIteration()
    else:
      it = m.iterations[n]
    iterations.append(it)
  for i,it in enumerate(iterations):
    ret += """
echo '**************************************************************'
echo '%s minimization'
echo '**************************************************************'
""" % ordinals[i]
    itparams = itparams0    
    if it.rcut != 1500 and len(grid_used) == 0: itparams += " --rcut %s" % str(it.rcut)
    if it.vmax != 100: 
      if it.mc:
        itparams += " --mcmax %s" % str(it.vmax)
      else:
        itparams += " --vmax %s" % str(it.vmax)
    if it.traj: itparams += " --traj"
    if it.prep: itparams += " --only-rot"
    if it.mc: 
      itparams += " --mc"      
      if it.mctemp != 3.5: itparams += " --mctemp %s" % str(it.mctemp)
      if it.mcscalerot != 0.05: itparams += " --mcscalerot %s" % str(it.mcscalerot)
      if it.mcscalecenter != 0.1: itparams += " --mcscalecenter %s" % str(it.mcscalecenter)
      if it.mcscalemode != 0.1: itparams += " --mcscalemode %s" % str(it.mcscalemode)
      if it.mcensprob != 0.1: itparams += " --mcensprob %s" % str(it.mcensprob)
    if has_restraints and it.restweight != 1: 
      itparams += " --restweight %s" % it.restweight
    
    if m.np > 1:
      attract = "python $ATTRACTDIR/../protocols/attract.py"
      tail = "$parals --output"  
    else:
      attract = "$ATTRACTDIR/attract"
      tail = ">"  
    if i == len(iterations) - 1:
      outp = "out_$name.dat"
    else:  
      outp = "stage%d_$name.dat" % (i+1)  
    outp0 = outp
    ipar="$params"
    if it.prep: 
      outp = outp0 + "-fixre"
      ipar="$paramsprep"
    gridpar = ""
    if len(gridparams) and not it.prep: gridpar = " $gridparams" 
    ret += "%s %s %s%s%s %s %s\n" % (attract, inp, ipar, gridpar, itparams, tail, outp0)
    if it.prep:
      ret += "$ATTRACTDIR/fix_receptor %s %d%s | python $ATTRACTTOOLS/fill.py /dev/stdin %s > %s\n" % (outp0, len(m.partners), flexpar1, outp0, outp)
    inp = outp       
  ret += "\n"  
  result = outp  
        
  if m.rescoring:
    ret += """
echo '**************************************************************'
echo 'Final rescoring'
echo '**************************************************************'
"""
    if m.np > 1:
      attract = "python $ATTRACTDIR/../protocols/attract.py"
      tail = "$parals --output"  
    else:
      attract = "$ATTRACTDIR/attract"
      tail = ">"  
    rcutsc = "--rcut %s" % str(m.rcut_rescoring)
    ret += "%s %s $scoreparams %s %s out_$name.score\n" \
     % (attract, result, rcutsc, tail)
    ret += """     
echo '**************************************************************'
echo 'Merge the scores with the structures'
echo '**************************************************************'
"""
    ret += "python $ATTRACTTOOLS/fill-energies.py %s out_$name.score > out_$name-scored.dat\n" % (result)
    ret += "\n"  
    result = "out_$name-scored.dat"

  if m.sort:
    ret += """
echo '**************************************************************'
echo 'Sort structures'
echo '**************************************************************'
"""
    ret += "python $ATTRACTTOOLS/sort.py %s > out_$name-sorted.dat\n" % result
    ret += "\n"  
    result = "out_$name-sorted.dat"

  if m.deredundant:    
    
    if m.fix_receptor == False or m.search != "syst": 
      ret += """
echo '**************************************************************'
echo 'Fix all structures into the receptor's coordinate frame
echo '**************************************************************'
"""     
      fixresult = result + "-fixre"
      ret += "$ATTRACTDIR/fix_receptor %s %d%s | python $ATTRACTTOOLS/fill.py /dev/stdin %s > %s\n" % (result, len(m.partners), flexpar1, result, fixresult)
      result = fixresult
    
    outp = os.path.splitext(result)[0] +"-dr.dat"
    ret += """
echo '**************************************************************'
echo 'Remove redundant structures'
echo '**************************************************************'
"""     
    ret += "$ATTRACTDIR/deredundant %s %d%s | python $ATTRACTTOOLS/fill-deredundant.py /dev/stdin %s > %s\n" % (result, len(m.partners), flexpar1, result, outp)
    ret += "\n"  
    result = outp

  if m.iattract is not None:
    ret += """
echo '**************************************************************'
echo 'iATTRACT refinement'
echo '**************************************************************'
"""                
    iattract_input = "out_$name-pre-iattract.dat"
    iattract_output = "out_$name-iattract.dat"
    iattract_output_demode = "out_$name-iattract-demode.dat"
    ret += "$ATTRACTTOOLS/top %s %d > %s\n" % (result, m.iattract.nstruc, iattract_input)
    iattract_params = ""
    iattract_params += " --infinite" #for now, use attract-infinite with iATTRACT
    iattract_params += " " + iattract_input
    iattract_params += " " + parameterfiledict["OPLSX"]
    for f in aa_filenames:
      iattract_params += " " + f            
    iattract_params += " --cdie --epsilon 10 --fix-receptor"
    iattract_params += " --icut %s" % m.iattract.icut
    iattract_params += " --np %s" % m.np
    iattract_params_demode = ""
    if ens_any:
      for pnr,p in enumerate(m.partners):
        f = aa_ensemble_lists[pnr]
        if f is not None:
          iattract_params += " --ens %d %s" % (pnr+1, f)
          iattract_params_demode += " --ens %d" % (pnr+1)
    if modes_any:
      iattract_params += " --modes %s" % aa_modesfile    
    if m.runname == "attract":
      iname = "i$name"
    else:
      iname = "iattract-$name"
    ret += "python $ATTRACTDIR/../protocols/iattract.py %s --name %s --output %s\n" % (iattract_params, iname, iattract_output)
    result = iattract_output
    if m.demode:
      ret += "python $ATTRACTTOOLS/demode.py %s %s > %s\n" % (iattract_output, iattract_params_demode, iattract_output_demode)    
      result = iattract_output_demode
    if m.fix_receptor or m.deredundant or m.calc_lrmsd:
      iattract_output_fixre = result + "-fixre"
      ret += "$ATTRACTDIR/fix_receptor %s %d%s | python $ATTRACTTOOLS/fill.py /dev/stdin %s > %s\n" % (result, len(m.partners), flexpar_iattract, result, iattract_output_fixre)        
      if m.fix_receptor or m.deredundant:
         result = iattract_output_fixre
  ret += """
echo '**************************************************************'
echo 'Soft-link the final results'
echo '**************************************************************'
"""       
  ret += "ln -s %s result.dat\n" % result

  result0 = result
  if m.collect:
    nr = m.nr_collect
    ret += """
echo '**************************************************************'
echo 'collect top %d structures:'
echo '**************************************************************'
""" % nr
    ret += "$ATTRACTTOOLS/top %s %d > out_$name-top%d.dat\n" % (result, nr, nr)
    collect_filenames = " ".join(aa_filenames)
    demodestr = ""
    flexpar_collect = ""
    if modes_any:      
      if not m.demode: 
        flexpar_collect = " --modes %s" % aa_modesfile 	  
      elif not m.iattract:        
        demodestr = "-demode"
        ret += "python $ATTRACTTOOLS/demode.py out_$name-top%d.dat > out_$name-top%d.dat-demode\n" % (nr, nr)        
    if m.iattract and not m.demode:
	flexpar_collect += " --name %s" %iname
    for pnr,p in enumerate(m.partners):
      if aa_ensemble_lists[pnr] is not None:
        flexpar_collect += " --ens %d %s" % (pnr+1,aa_ensemble_lists[pnr])    
    ret += "$ATTRACTDIR/collect out_$name-top%d.dat%s %s%s > out_$name-top%d.pdb\n" % \
     (nr, demodestr, collect_filenames, flexpar_collect, nr)
    ret += "ln -s out_$name-top%d.pdb result.pdb\n" % nr
    ret += "\n"
      
  if m.deflex:
    deflex_header = """
echo '**************************************************************'
echo 'Remove flexibility for RMSD calculations'
echo '**************************************************************'
tmpf=`mktemp`
tmpf2=`mktemp`

""" 
    outp, outp2 = "$tmpf", "$tmpf2"
    if modes_any and not (m.iattract and m.demode):           
      ret += deflex_header
      has_tmpf = True
      ret += "python $ATTRACTTOOLS/demode.py %s > %s\n" % \
        (result, outp)
      result = outp
      outp, outp2 = outp2, outp          
    if ens_any:
      if not has_tmpf:
        ret += deflex_header
        has_tmpf = True
      ret += "python $ATTRACTTOOLS/de-ensemblize.py %s > %s\n" % \
        (result,outp)
      result = outp
      outp, outp2 = outp2, outp

  if m.calc_lrmsd or m.calc_irmsd or m.calc_fnat:

    flexpar2 = ""        
    if not m.deflex and not m.demode and m.iattract is not None:
      flexpar2 += " --name %s" % iname
    if modes_any and not m.demode and not m.deflex: 
      flexpar2 = " --modes %s" % aa_modesfile    
    for pnr,p in enumerate(m.partners):
      if m.deflex == False and aa_ensemble_lists[pnr] is not None:
        flexpar2 += " --ens %d %s" % (pnr+1,aa_ensemble_lists[pnr])
    
    rmsd_refenames = [None] * len(m.partners)
    for pnr in range(len(m.partners)):
      filename = rmsd_filenames[pnr]
      p = m.partners[pnr]        
      if p.rmsd_pdb is not None:
        filename = p.rmsd_pdb.name
      rmsd_refenames[pnr] = filename
        
    bb_str = "backbone"
    if m.rmsd_atoms == "all":
      bb_str = "all-atom"
    elif m.rmsd_atoms == "trace":
      mt = m.partners[0].moleculetype
      if mt == "Protein": bb_str = "c-alpha"
      elif mt in ("DNA", "RNA"): bb_str = "phosphate"
        
  if m.calc_lrmsd:      
    ret += """
echo '**************************************************************'
echo 'calculate %s ligand RMSD'
echo '**************************************************************'
""" % bb_str      
  
    if m.fix_receptor == False and (not m.deredundant or m.iattract): 
      fixresult = result + "-fixre"
      ret += "$ATTRACTDIR/fix_receptor %s %d%s | python $ATTRACTTOOLS/fill.py /dev/stdin %s > %s\n" % (result, len(m.partners), flexpar2, result, fixresult)
      result = fixresult
  
    lrmsd_allfilenames = []
    for f1, f2 in zip(rmsd_filenames[1:], rmsd_refenames[1:]):
      lrmsd_allfilenames.append(f1)
      lrmsd_allfilenames.append(f2)
    lrmsd_allfilenames = " ".join(lrmsd_allfilenames)
    lrmsdpar = "--receptor %s" % rmsd_filenames[0]
    if m.rmsd_atoms == "all":
      lrmsdpar += " --allatoms"
    elif m.rmsd_atoms == "trace":
      mt = m.partners[0].moleculetype
      if mt == "Protein": lrmsdpar += " --ca"
      elif mt in ("DNA", "RNA"): lrmsdpar += " --p"
    lrmsdresult = os.path.splitext(result0)[0] + ".lrmsd"
    ret += "python $ATTRACTDIR/lrmsd.py %s %s%s %s > %s\n" % (result, lrmsd_allfilenames, flexpar2, lrmsdpar, lrmsdresult)
    ret += "ln -s %s result.lrmsd\n" % lrmsdresult
    ret += "\n"

  if m.calc_irmsd:      
    ret += """
echo '**************************************************************'
echo 'calculate %s interface RMSD'
echo '**************************************************************'
""" % bb_str      
    
    irmsd_allfilenames = []
    for f1, f2 in zip(rmsd_filenames, rmsd_refenames):
      irmsd_allfilenames.append(f1)
      irmsd_allfilenames.append(f2)
    irmsd_allfilenames = " ".join(irmsd_allfilenames)
    irmsdresult = os.path.splitext(result0)[0] + ".irmsd"
    bbo = "" 
    if m.rmsd_atoms == "all": bbo = "--allatoms"
    ret += "python $ATTRACTDIR/irmsd.py %s %s%s %s > %s\n" % (result, irmsd_allfilenames, flexpar2, bbo, irmsdresult)
    ret += "ln -s %s result.irmsd\n" % irmsdresult
    ret += "\n"

  if m.calc_fnat:      
    ret += """
echo '**************************************************************'
echo 'calculate fraction of native contacts'
echo '**************************************************************'
"""    
    fnat_allfilenames = []
    for f1, f2 in zip(rmsd_filenames, rmsd_refenames):
      fnat_allfilenames.append(f1)
      fnat_allfilenames.append(f2)
    fnat_allfilenames = " ".join(fnat_allfilenames)
    fnatresult = os.path.splitext(result0)[0] + ".fnat"
    ret += "python $ATTRACTDIR/fnat.py %s 5 %s%s > %s\n" % (result, fnat_allfilenames, flexpar2, fnatresult)
    ret += "ln -s %s result.fnat\n" % fnatresult
    ret += "\n"

  if m.calc_lrmsd or m.calc_irmsd or m.calc_fnat:
    if has_tmpf:
      ret += "rm -f $tmpf $tmpf2\n"
      result = result0
    if m.calc_lrmsd:
      if m.fix_receptor == False:
        ret += "rm -f %s\n" % fixresult
    
  if len(cleanupfiles):  
    ret += "\nrm -f " + " ".join(cleanupfiles) + "\n"
  if (m.footer is not None):
    ret += m.footer + "\n\n"
    
  ret += "fi ### move to disable parts of the protocol\n"
  ret = ret.replace("\n\n\n","\n\n")
  ret = ret.replace("\n\n\n","\n\n")
  ret = ret.replace("\r\n", "\n")
  return ret

