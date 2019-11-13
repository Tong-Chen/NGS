#! /usr/bin/env python

import pygame
import math

#   version 1.2.2

# Arbitarary list of gel sizes for easy reference
# Tuples give value of height in pixels and number of lanes
GEL_SIZE = {'small': (240, 4),
            'big': (360, 8),
            'long single': (300, 1)}

# Approximate values for 1kb ladder
# Note that samples are list of tuples in the form (length , concentration)
LADDER_1KB =[
    (300, 250),
    (500, 150),
    (700, 100),
    (1000, 75),
    (1500, 50),
    (2000, 40),
    (2500, 35),
    (3000, 30),
    (4000, 25),
    (5000, 30),
    (6000, 15),
    (8000, 12),
    (10000, 15),
    (14000, 4),
    (24000, 2)]

def gamma(x, t, k=1):
    b = 1.0 / t
    x = x * k
    c = 0.0

    for i in range(0, k):
        c += (math.exp(-b * x) * (b * x) ** i) / math.factorial(i)

    return c

def testColour(c):
    """ Convert integer to colour, which saturates if > 255 """

    if c <= 255:
        return (c, c, c)
    else:
        return (255,255,255)

def display(width, height, image):
    """ Create a Pygame image of dimensions width x height and show image """

    screen = pygame.display.set_mode((width, height))
    screen.blit(image, (60, 30))
    pygame.display.flip()

    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

class Gel:
    """ A Gel contains lanes into which DNA can be loaded and run. 
        The Gel is exposured to see the location and intensity of DNA. """

    LANE_WIDTH = 30.
    LANE_HEIGHT = 4
    LANE_MARGIN = 6
    BORDER = 12

    def __init__(self, size, agarose=1):
        self.y = GEL_SIZE[size][0]
        self.x = 2 * Gel.BORDER + GEL_SIZE[size][1] * (Gel.LANE_WIDTH + Gel.LANE_MARGIN) - Gel.LANE_MARGIN

        self.agarose = agarose      # Should test whether value is reasonable (say, 0.5% - 3%)
        self.optimum_DNA_length = 2000 / agarose ** 3   # Best separate based on agarose hole size

        self.samples = []
        self.lanes = []

        for n in range(GEL_SIZE[size][1]):
            self.samples.append([])
            self.lanes.append([])

        for lane in self.lanes:
            for n in range(Gel.BORDER, self.y + 3):
                lane.append(0.0)
        print "*****Initialzie****"
        print self.samples
        print self.lanes

    def loadSample(self, lane, sample):
        """Add list containing tuple of DNA length and concentrations to lane in gel. """

        for dna in sample:
            strand = {'length': float(dna[0]), 'conc': dna[1] * 1./Gel.LANE_HEIGHT, 'position': Gel.BORDER+1}
            self.samples[lane-1].append(strand)
        print "***loadSample***",
        print lane,
        print sample
        print self.samples
        print self.lanes

    def run(self, time=30.0, voltage=20.0):
        """ Move loaded DNA down the gel at a rate dependent on voltage, DNA length and agarose concentration. """

        max_dist = 0.25 * time * voltage
        for sample in self.samples:
            for dna in sample:
                g = gamma( dna['length']/20, int(self.optimum_DNA_length/20) )
                dna['position'] += max_dist * g
        print self.samples
        print self.lanes

    def expose(self, exposure=0.1, aperture=2):
        """Returns an image of the gel with the DNA highlighted"""

        c1, c2 = exposure * 100, exposure * 200
        fill_colour = (c1, c1, c2)
        image = pygame.Surface((self.x+2, self.y+2))
        pygame.draw.rect(image, fill_colour, (1,1, self.x, self.y), 0)

        edge_colour = (140, 140, 150)
        edges = pygame.Surface((self.x+2, self.y+2))
        edges.set_colorkey((0,0,0))
        edges.set_alpha(120)
        pygame.draw.rect(image, edge_colour, (0,0, self.x+2, self.y+2), 1)

        self._findDNAConcentrations(c1)

        x = Gel.BORDER + 1

        for lane in self.lanes:
            for position in range(len(lane)):
                if lane[position] > 0:

                    brightness = lane[position] * exposure / aperture + c1
                    colour1 = testColour(brightness)
                    colour2 = testColour(brightness*0.6)
                    colour3 = testColour(brightness*0.3)

                    pygame.draw.line(image, colour3, (x-1, position+Gel.BORDER),(x+Gel.LANE_WIDTH+1, position+Gel.BORDER))
                    pygame.draw.line(image, colour2, (x, position+Gel.BORDER),(x+Gel.LANE_WIDTH, position+Gel.BORDER))
                    pygame.draw.line(image, colour1, (x+1, position+Gel.BORDER),(x+Gel.LANE_WIDTH-1, position+Gel.BORDER))
    
            pygame.draw.rect(edges, edge_colour, (x, Gel.BORDER, Gel.LANE_WIDTH, Gel.LANE_HEIGHT), 1)
            x += Gel.LANE_WIDTH + Gel.LANE_MARGIN

        image.blit(edges, (0,0))
        return image

    def _findDNAConcentrations(self, background):
        """Determines where in the concentration of DNA in every part of the gel"""

        length = len(self.lanes[0])

        for x in range(len(self.samples)):
            for dna in self.samples[x]:
                for y in range(Gel.LANE_HEIGHT-2):
                    pos = int(dna['position']) + y
                    if pos < length-4:
                        # Very crude way to create blurred line
                        self.lanes[x][pos-2] += 0.06 * dna['conc'] * dna['length']
                        self.lanes[x][pos-1] += 0.12 * dna['conc'] * dna['length']
                        self.lanes[x][pos] += 0.2 * dna['conc'] * dna['length']
                        self.lanes[x][pos+1] += 0.12 * dna['conc'] * dna['length']
                        self.lanes[x][pos+2] += 0.06 * dna['conc'] * dna['length']
                    print self.lanes

def example():
    """ An example of how to set up and run a gel. """

    # Samples are list of tuples in the form (length, concentration)
    sample1 = [(400, 200), (500, 200), (900, 200)]
    sample2 = [(400, 200), (900, 200), (1500, 200)]
    sample3 = [(4000, 20), (4100, 20), (4900, 10)]

    # Set up gel, giving its size and % agarose
    myGel = Gel(size="small", agarose=1.2)

    # Load samples
    myGel.loadSample(1, LADDER_1KB)
    myGel.loadSample(2, sample1)
    myGel.loadSample(3, sample2)
    myGel.loadSample(4, sample3)     
    myGel.loadSample(4, sample1) # you can load more than 1 sample per lane

    # Run gel for 45 minutes
    myGel.run(time=45)

    # Take a picture of your gel with a 50ms exposure
    my_image = myGel.expose(exposure=0.05)

    # Display image on pygame screen
    #display(320, 320, my_image)

    # Or save image as a PNG
    pygame.image.save(my_image, 'test_run.png')

def main():
    example()

if __name__ == '__main__':
    main()
