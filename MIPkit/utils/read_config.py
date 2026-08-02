import os
import yaml 

def read_config(args):
    # ROOT_DIR will be utils, go up one
    ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    conf = f"{ROOT_DIR}/config/config.yaml"
    if os.path.isfile(conf):
        with open(conf) as stream:
            yamlin = yaml.safe_load(stream)

            if yamlin.get("obabel"):
                args.obabel = yamlin["obabel"]

            if "docking-method" in yamlin:
                args.dockmethod = yamlin["docking-method"]

            if yamlin.get("vina"):
                args.vina = yamlin["vina"]

            if yamlin.get("gnina"):
                args.gnina = yamlin["gnina"]

            if yamlin.get("acpype"):
                args.acpype = yamlin["acpype"]

            if yamlin.get("gmx"):
                args.gmx = yamlin["gmx"]
        return args

    else:
        print(
            f"No config file {conf} available. We need one to direct the scripts properly. Exiting... "
        )
        quit()