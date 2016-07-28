from orangecontrib.bio import resolwe

gen = resolwe.connect('anonymous@genialis.com', 'anonymous', 'https://dictyexpress.research.bcm.edu', 'genesis')

experiments = gen.fetch_etc_objects()

for exp in experiments:
    print("Experiment id: " + exp.id + " - " + str(exp))
