import math

LAMINAR_BOUNDARY = 2100
ACCELERATION_DUE_TO_GRAVITY = 9.80665


class Segment:
    
    def __init__(self, length=0, K_fittings=0, elevation_gain=0,
                 mass_flowrate=0, fluid=None, pipe=None, node_in = None,
                 node_out = None):
        self.length = length
        self.K_fittings = K_fittings # TODO replace with elements dict
        self.elevation_gain = elevation_gain
        self.mass_flowrate = mass_flowrate
        self.fluid = fluid
        self.pipe = pipe
        self.node_in = node_in
        self.node_out = node_out
        
        self.cross_sectional_area = math.nan
        self.volumetric_flowrate = math.nan
        self.velocity = math.nan
        self.reynolds_number = math.nan
        self.friction_factor = math.nan
        self.K_pipe = math.nan
        self.K_total = math.nan
        self.pressure_drop_static = math.nan
        self.pressure_drop_pipe = math.nan
        self.pressure_drop_elements = math.nan
        self.pressure_drop_dynamic = math.nan
        self.pressure_drop_total = math.nan
    
        #self.update()
        
        
    def update(self):
        """Calculate pipe pressure drop."""
        
        seg = (self.node_in, self.node_out)
        
        if self.fluid is None:
            print("Warning: Fluid not specified for segment {}".format(seg))
        if self.pipe is None:
            print("Warning: Pipe not specified for segment {}".format(seg))
            
        self.cross_sectional_area = (0.25 * math.pi * self.pipe.diameter *
                                     self.pipe.diameter)
        
        if self.fluid.density > 0:
            self.volumetric_flowrate = self.mass_flowrate / self.fluid.density
        else:
            self.volumetric_flowrate = math.nan
        
        if self.cross_sectional_area > 0:
            self.velocity = (self.volumetric_flowrate /
                             self.cross_sectional_area)
        else:
            self.velocity = math.nan
        
        if self.fluid.viscosity > 0:
            self.reynolds_number = (self.fluid.density * self.velocity *
                                    self.pipe.diameter / self.fluid.viscosity)
            if self.reynolds_number < LAMINAR_BOUNDARY:
                print("Warning: Laminar flow in segment {}".format(seg))
        else:
            self.reynolds_number = math.nan
        
        self.friction_factor = friction_factor(self.reynolds_number,
                                               self.pipe.relative_roughness)
        
        self.pressure_drop_static = (self.fluid.density * self.elevation_gain *
                                     ACCELERATION_DUE_TO_GRAVITY)
        
        if self.diameter > 0:
            self.K_pipe = (self.friction_factor * self.length /
                           self.pipe.diameter)
        else:
            self.K_pipe = math.nan
            
        self.K_total = self.K_pipe + self.K_fittings
        
        self.pressure_drop_pipe = (self.K_pipe * 0.5 * self.fluid.density *
                                   self.velocity * self.velocity)


        self.pressure_drop_elements = (self.K_fittings * 0.5 *
                                       self.fluid.density *
                                       self.velocity * self.velocity)
        
        self.pressure_drop_dynamic = (self.pressure_drop_pipe +
                                      self.pressure_drop_elements)
        
        self.pressure_drop_total = (self.pressure_drop_static +
                                    self.pressure_drop_dynamic)


class Fluid:
    def __init__(self, density, viscosity):
        self.density = density
        self.viscosity = viscosity


class Pipe:
    def __init__(self, diameter, roughness, nominal_diameter="",
                 pipe_class="", schedule=""):
        self.diameter = diameter
        self.roughness = roughness
        self.nominal_diameter = str(nominal_diameter)
        self.pipe_class = str(pipe_class)
        self.schedule = str(schedule)
        self.relative_roughness = self.roughness / self.diameter


def friction_factor(reynolds_number, relative_roughness):
    """Calculate the Darcy friction factor."""
    # Poiseuille's Law for laminar flows.
    if reynolds_number < LAMINAR_BOUNDARY: return reynolds_number / 64
    
    # Colebrook-White equaiton for turbulent flows.
    # Goudar-Sonnad correlation used for explicit calculation.
    a = 2 / math.log(10)
    b = relative_roughness / 3.7
    d = math.log(10) * reynolds_number / 5.02
    s = b * d + math.log(d)
    q = math.pow(s, (s / (s + 1)))
    g = b * d + math.log(d / q)
    z = math.log(q / g)
    D_LA = z * g / (g + 1)
    D_CFA = D_LA * (1 + (z / 2) / (math.pow(g + 1,2) + (z / 3) * (2 * g - 1)))
    return math.pow(a * (math.log(d / q) + D_CFA), -2)


if __name__ == "__main__":
    
    fluid = Fluid(998, 0.001)
    pipe = Pipe(0.05, 0.000045, 50)
    
    segment = Segment()
    segment.fluid = fluid
    segment.pipe = pipe
    segment.length = 10
    segment.diameter = 0.05
    segment.roughness = 0.000045
    segment.mass_flowrate = 4
    
    segment.update()
    
    variables = vars(segment)
    for v in variables:
        print('{0:26}:  {1}'.format(v, variables[v]))
