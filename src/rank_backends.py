import subprocess
import os
import os.path
import tempfile

class RankBackend(object):
    pass

class MatlabBackend(RankBackend):
    def compute_ranks(self,input_dir,output_file):
        matlabenv = dict(os.environ);
        matlabenv["MATLABPATH"]=self._matlab_path
        try:
            with open(output_file,'w') as outf:
                for infile in os.listdir(input_dir):
                    base_name = infile.split('.')[0]
                    parts = base_name.split('_')
                    md = tuple(map(int,parts[1:]))
                    with tempfile.NamedTemporaryFile() as temp:
                        with open('/dev/null','w') as nullfile:
                            subprocess.check_call(["matlab","-nodesktop","-nosplash","-nodisplay","-r","{}('{}','{}');exit;".format(self._matlab_func,os.path.join(input_dir,infile),temp.name)],stdout=nullfile,env=matlabenv)
                        temp.seek(0)
                        rank = int(temp.readline().split()[0])
                        outf.write("{} {}\n".format(md,rank))
        except:
            os.remove(output_file)
            raise

#TODO write the specific versions of MatlabBackend

class MatlabQRBackend(MatlabBackend):
    def __init__(self):
        self._matlab_path = "./src/"
        self._matlab_func = "qr_single"

class SageBackend(RankBackend):
    def compute_ranks(self,input_dir,output_file):
        with open(output_file,'w') as output:
            subprocess.check_call(["sage","-q","./src/ranks.sage",input_dir],
                                  stdout=output)
