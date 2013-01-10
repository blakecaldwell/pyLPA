from matplotlib.colors import LinearSegmentedColormap

# yellow-red-black-blue-cyan colormap
yrbbcdict = {'red' : [(0.0, 1.0, 1.0),
                    (0.25, 1.0, 1.0),
                    (0.5, 0.0, 0.0),
                    (1.0, 0.0, 0.0)],
           'green' : [(0.0, 1.0, 1.0),
                      (0.25, 0.0, 0.0),
                      (0.75, 0.0, 0.0),
                      (1.0, 1.0, 1.0)],
           'blue': [(0.0, 0.0, 0.0),
                    (0.5, 0.0, 0.0),
                    (0.75, 1.0, 1.0),
                    (1.0, 1.0, 1.0)]
           }
boston = LinearSegmentedColormap('my_yrbbc', yrbbcdict, 256)

#reversed version
yrbbcdict_r = {'red' : [(0.0, 0.0, 0.0),
                    (0.5, 0.0, 0.0),
                    (0.75, 1.0, 1.0),
                    (1.0, 1.0, 1.0)],
           'green' : [(0.0, 1.0, 1.0),
                      (0.25, 0.0, 0.0),
                      (0.75, 0.0, 0.0),
                      (1.0, 1.0, 1.0)],
           'blue': [(0.0, 1.0, 1.0),
                    (0.25, 1.0, 1.0),
                    (0.5, 0.0, 0.0),
                    (1.0, 0.0, 0.0)]
           }

boston_r = LinearSegmentedColormap('my_yrbbc_r', yrbbcdict_r, 256)
