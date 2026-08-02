import yaml
import os

class init_attackers:
	def __init__(self):
		self.dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../initiators/"))

	def load(self, name):
		
		with open(f"{self.dir}/initiators.yaml") as stream:
		    try:
		        inits = yaml.safe_load(stream)

		    except yaml.YAMLError as exc:
		        print(exc)
		        exit()

		initlist = inits[name]

		atoms = []
		nums = []

		if initlist is not None:
			for i, atm in enumerate(initlist):
				if i % 2:
					nums.append(atm)
				else:
					atoms.append(atm)

		return atoms, nums
	
	def all(self):
		with open(f"{self.dir}/initiators.yaml") as stream:
		    try:
		        inits = yaml.safe_load(stream)

		    except yaml.YAMLError as exc:
		        print(exc)
		        exit()

		key = []
		for lbl in inits:
			key.append(lbl)

		return key


# init_attackers()
# print(init_attackers.load('APS'))