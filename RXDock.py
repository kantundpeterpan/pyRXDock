#!/usr/bin/env python
# coding: utf-8

import os
from multiprocessing import Pool
from subprocess import run
from glob import glob
import shlex
import re
from functools import partial


# In[52]:

class RXDock():

    @classmethod
    def rbdock(cls, ligand_file, receptor_prm, 
                  output_suffix = '_out', dock_prm = 'dock.prm',
                  **kwargs):
        """

        Parameters
        ----------
        ligand_file : str
            Path to SDF file containing ligands to dock
        receptor_prm : str
            Path to receptor Paramters
        output_suffix : str, optional
            Resulting sdf files will be append with this suffix.
            The default is '_out'.
        dock_prm : str, optional
            Name of docking protocol parameter file.
            The default is 'dock.prm'.
        **kwargs : arg:value
            all commandline flags for rbdock can be passed as
            keyword arguments. A simple flag is passed with an
            empty string as value, e.g. H='' (will append the -H flag)

        Returns
        -------
        str
            Absolute path to the sdf file with docked ligands.

        """
    
        run(
            shlex.split(
                    'rbdock -i %s -o %s -r %s -p %s %s' %
                    (
                        ligand_file,
                        ligand_file.replace('.sd', '')+output_suffix,
                        receptor_prm,
                        dock_prm,
                        ' '.join(['-%s %s' % (k,v) for k,v in kwargs.items()])
                        )
                )
        )
        
        return os.path.abspath(
            os.path.join(
                    os.path.dirname(ligand_file),
                        os.path.splitext(os.path.basename(ligand_file))[0])+\
                            output_suffix+'.sd')
    
    
    # In[53]:
    @classmethod
    def sdsplit(cls, ligand_file, dir_prefix = 'tmp',
                   split_prefix = 'tmp_',
                   n = 30):
        """
        Wraps the sdsplit tool

        Parameters
        ----------
        ligand_file : str
            path to SDF file to split.
        dir_prefix : str, optional
            split SDF files will be put into ./{dir_prefix}.
            The default is 'tmp'.
        split_prefix : str, optional
            basename for splitted SDF files.
            The default is 'tmp_'.
        n : int, optional
            Number of ligands in split files.
            The default is 30.

        Returns
        -------
        splits : list[str]
            Paths to splitted SDF files (./{dir_prefix}/{split_prefix}*.sd)

        """        
        if not os.path.exists(dir_prefix):
            os.mkdir(dir_prefix)
            
        run(
                [
                    'sdsplit',
                    '-%s' % str(n),
                    '-o %s/%s' % (dir_prefix, split_prefix),
                    '%s' % ligand_file
                    ]
            )
    
        temp = glob('%s/%s*.sd' % (dir_prefix, split_prefix))
    
        r = re.compile('\d+.sd')
        splits = sorted(temp, key = lambda x: int(r.findall(x)[0].replace('.sd','')))
    
        return splits   
    
    @classmethod
    def _multidock(cls, ligands, receptor_prm, output_suffix = '_out', 
               dock_prm = 'dock.prm', n_jobs = 2, **kwargs):
    
        p = Pool(n_jobs)
        
        _map_func = partial(cls.rbdock, receptor_prm = receptor_prm,
                            output_suffix=output_suffix, dock_prm=dock_prm,
                            **kwargs)
        
        out = p.map(_map_func, ligands)
        
        return out
    
    @classmethod
    def splitdock(cls, ligands, receptor_prm,
                  tmpdir_prefix = 'tmp',
                  split_prefix = 'tmp_',
                  n_splits = 30,
                  output_suffix = '_out',
                  dock_prm = 'dock.prm',
                  **kwargs):
        
        splitted = cls.sdsplit(ligands,
                               dir_prefix=tmpdir_prefix,
                               split_prefix=split_prefix,
                               n = n_splits)
        
        res = cls._multidock(splitted,
                             receptor_prm=receptor_prm,
                             output_suffix=output_suffix,
                             dock_prm=dock_prm,
                             n_jobs = kwargs.pop('n_jobs', 2),
                             **kwargs)
        
        return res
        