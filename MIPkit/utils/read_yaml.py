import yaml
import os

def read_yaml(args, recipe_only=False):
    if os.path.isfile(args.config[0]):
        with open(args.config[0]) as stream:
            yamlin = yaml.safe_load(stream)

            if yamlin.get("fms"):
                args.fms = yamlin["fms"]

            if yamlin.get("protein") and not recipe_only:
                args.protein = [yamlin["protein"]]

            if yamlin.get("scale") and not recipe_only:
                args.scale = yamlin["scale"]

            if yamlin.get("center") and not recipe_only:
                args.center = yamlin["center"]

            if yamlin.get("size") and not recipe_only:
                args.size = yamlin["size"]

            if yamlin.get("energy_range") and not recipe_only:
                args.energy_range = yamlin["energy_range"]

            if yamlin.get("exhaustiveness") and not recipe_only:
                args.exhaustiveness = yamlin["exhaustiveness"]

        return args

    else:
        print(f"No config file {args.config[0]} available. Exiting... ")
        quit()