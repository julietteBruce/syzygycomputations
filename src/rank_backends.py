import subprocess

class RankBackend(object):
    pass

class MatlabBackend(RankBackend):
    def compute_ranks(self,input_dir,output_file):
        subprocess.check_call(["matlab","-nodesktop","-nosplash","-nodisplay","-r {}({},{})".format(self.__matlab_func,input_dir,output_file)])

#TODO write the specific versions of MatlabBackend
                               
class SageBackend(RankBackend):
    def compute_ranks(self,input_dir,output_file):
        with open(output_file,'w') as output:
            subprocess.check_call(["sage","-q","./src/ranks.sage",input_dir],
                                  stdout=output)
