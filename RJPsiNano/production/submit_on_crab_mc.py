from CRABClient.UserUtilities import config, ClientException
#from input_crab_data import dataset_files
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser
import os 
production_tag = datetime.date.today().strftime('%Y%b%d')


config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'RJPsiNANO_%s' % production_tag

config.section_('Data')
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/cgalloni/%s' % ('crab_job_' + production_tag)
config.Data.inputDBS = 'global'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/run_nano_jpsi_2mu_cfg.py'
config.JobType.inputFiles = ['crab_script_mc_job.sh', '../test/run_nano_jpsi_2mu_cfg.py', 'puReweight_2016.py', 'pileup/pileup_2016GH.root']
#'pileup/PileupHistogram-UL2016-100bins_withVar.root']
#config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script_mc_job.sh'
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CSCS'

if __name__ == '__main__':

  os.system("cp crab_script_mc.sh crab_script_mc_job.sh; sed -i s/TAGDATE/{}/ crab_script_mc_job.sh;".format(production_tag))
  from CRABAPI.RawCommand import crabCommand
  from CRABClient.ClientExceptions import ClientException
  from httplib import HTTPException
  from multiprocessing import Process

  def submit(config):
      try:
          crabCommand('submit', config = config)
      except HTTPException as hte:
          print "Failed submitting task: %s" % (hte.headers)
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)


  parser = ArgumentParser()
  parser.add_argument('-y', '--yaml', default = 'samples_mc_rjpsi_2016.yml', help = 'File with dataset descriptions')
  parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed')
  args = parser.parse_args()

  with open(args.yaml) as f:
    doc = yaml.load(f) # Parse YAML file
    common = doc['common'] if 'common' in doc else {'data' : {}, 'mc' : {}}
    
    # loop over samples
    for sample, info in doc['samples'].iteritems():
      # Given we have repeated datasets check for different parts
      parts = info['parts'] if 'parts' in info else [None]
      for part in parts:
        name = sample % part if part is not None else sample
        
        # filter names according to what we need
        if not fnmatch(name, args.filter): continue
        print 'submitting', name

        isMC = info['isMC']
        config.Data.inputDataset = info['dataset'] % part \
                                   if part is not None else \
                                      info['dataset']

        #config.General.requestName = name + '2'
        config.General.requestName = name 
        common_branch = 'mc' if isMC else 'data'
        #config.Data.splitting = 'FileBased' if isMC else 'LumiBased'
        config.Data.splitting = 'EventAwareLumiBased' if isMC else 'LumiBased' 
        config.Data.unitsPerJob = 15000  #for 2017  
        if not isMC:
            config.Data.lumiMask = info.get(
                'lumimask', 
                common[common_branch].get('lumimask', None)
            )
        else:
            config.Data.lumiMask = ''

        #config.Data.unitsPerJob = info.get(
        #    'splitting',
        #    common[common_branch].get('splitting', None)
        #)
        globaltag = info.get(
            'globaltag',
            common[common_branch].get('globaltag', None)
        )
        
        config.JobType.pyCfgParams = [
            'isMC=%s' % isMC, 'reportEvery=5000',
            'tag=%s' % production_tag,
            'globalTag=%s' % globaltag,
        ]
        
        config.JobType.outputFiles = ['_'.join(['RJPsi', 'mc' if isMC else 'data','pu', production_tag])+'.root']
        #config.JobType.outputFiles = ['_'.join(['RJPsi', 'mc' if isMC else 'data', production_tag])+'.root']
        
        print config
        submit(config)

