from math import sqrt, pi
import pygame
from random import random
import time
import os
from PIL import Image
import numpy as np


class BarnsHut:
    def __init__(self, particles, a, b, c=2, Θ=0.0, color=(255, 0, 0)):
        """

        :param particles: particles contained
        :param a: first corner of boundry box
        :param b: second corner of boundry box
        :param c: divisions per dimension
        :param Θ: distance to size ratio
        """
        # if No particles: no mass, and center at average point
        self.color = color
        self.particles = particles
        self.Θ = Θ
        self.a = a
        self.b = b
        self.c = c
        self.M = 0.0
        self.s = np.power(b.delta(a).inner_mul(), 1.0/self.a.dim)
        self.location = (a + b) * 0.5
        self.children = []
        if len(particles) > 0:
            # if there are particles: find center of mass and total mass
            self.location = Vector([0, 0])
            self.M = 0.0
            momentum = Vector([0, 0])
            for p in particles:
                self.location.translate(p.location * p.M)
                self.M += p.M
            self.location.div(self.M)
            if len(particles) > 1:
                beg = time.time()
                self.children = self.create_children(random_color())
                #   print("time to create children: " + str(time.time()-beg))
            """else:
                draw_rect(self.a, self.b, self.color)"""
        self.R = sqrt(self.M/(pi*P))
        #   print("p: " + str(len(self.particles)) + ", c: " + str(len(self.children)))
        #   pygame.draw.circle(screen, (255, 0, 255), self.location.get_translate(relative_pos).to_int_list(), int(self.R))


    def create_children(self, color):
        children = []
        child_count = self.c**self.a.dim    # amount of children
        skip = self.b.delta(self.a).div(self.c)
        for num in range(child_count):
            s = self.a.copy()
            tmp = num
            loc = []
            # you can see the space as a number, in base [c], with [dimensions] digits.
            for d in range(s.dim):
                loc.append(tmp % self.c)
                s.add_to_axis(skip.coords[d]*(tmp % self.c), d)
                tmp = tmp // self.c
            #   print(loc)
            children.append(BarnsHut(self.particles_in_area(s, s+skip), s, s+skip, color=color))
        return children

    def particles_in_area(self, s, t):
        good = []
        for p in self.particles:
            if p.within_boundry(s, t):
                good.append(p)
        return good

    def enforce(self, particle):
        how_many_particles = len(self.particles)
        if how_many_particles == 1:
            if not(particle is self.particles[0]):
                if particle.feel_force_from(self.particles[0], can_merge=True):
                    self.particles = []
        elif how_many_particles > 1:
            if not self.distant(particle) or (particle in self.particles):
                for c in self.children:
                    c.enforce(particle)
            else:
                particle.feel_force_from(self)

    def distant(self, p):
        if self.R*2 / (p.distance(self) + 0.000000001) <= self.Θ:
            return True

    def delta(self, p2):
        return self.location.delta(p2.location)

    def depth(self, n):
        if len(self.children) > 0:
            return max([a.depth(n+1) for a in self.children])
        return n

    def nodes(self):
        return 1 + sum([a.nodes() for a in self.children])


class Vector:
    def __init__(self, coords):
        self.coords = np.array(coords).astype(float)
        self.dim = sum(self.coords.shape)

    def translate(self, p2):
        self.coords += p2.coords
        return self

    def inner_mul(self):
        return np.prod(self.coords)

    def within_boundry(self, a, b):
        return np.all(self.coords <= b.coords) and np.all(self.coords > a.coords)

    def add_to_axis(self, n, d):
        self.coords[d] += n

    def __add__(self, p2):
        return Vector(self.coords + p2.coords)

    def get_translate(self, p2):
        return Vector(self.coords + p2.coords)

    def multiply(self, n):
        self.coords *= n
        return self

    def __mul__(self, n):
        return Vector(self.coords * n)

    def get_multiply(self, n):
        return Vector(self.coords * n)

    def div(self, n):
        self.coords /= n
        return self

    def get_div(self, n):
        return Vector(self.coords / n)

    def copy(self):
        return Vector(np.copy(self.coords))

    def distance(self, p2):
        return np.sqrt(np.sum((self.coords-p2.coords)**2))

    def delta(self, p2):
        return Vector(self.coords - p2.coords)

    def to_int_list(self):
        return self.coords.astype(int)


class Force:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

    def calc_force(self):
        # F = GM1M2/R²
        R = self.p1.distance(self.p2)

        R_eff = max(R, self.p1.R + self.p2.R)
        size_of_F = G * self.p1.M * self.p2.M / (R_eff ** 2)
        direction = self.p2.delta(self.p1).div(R_eff)  # relative to p1. For p2, it's -1 times that. vector with direction of force, size=1
        return direction.multiply(size_of_F), R  # Force vector: size of force, direction of force relative to p1

    def enforce(self):
        for p in [self.p1, self.p2]:
            if p not in Particles:
                self.p1 = None
                self.p2 = None
                return
        F, R = self.calc_force()
        if F:
            self.p1.feel_force(F)
            F.multiply(-1.0)  # flip direction of force to make it relative to p2
            self.p2.feel_force(F)
        if R <= COLLISION_COEFFICIENT*(self.p1.R + self.p2.R):
            # create a new particle containing both masses, found at center of mass
            Particles.remove(self.p1)
            Particles.remove(self.p2)
            add_this_particle(total_mass([self.p1, self.p2]))

            self.p1 = None
            self.p2 = None
            return



class Particle:
    def __init__(self, M, P, color=(0, 0, 0), location=Vector([0, 0]), V=None, Moment=None):
        self.color = color
        self.R = sqrt(M/(pi*P))
        self.P = P
        self.M = M
        self.location = location
        if V:
            self.V = V
        elif Moment:
            self.V = Moment.div(self.M)
        else:
            self.V = Vector([0, 0])
        self.ΣF = Vector([0, 0])

    def within_boundry(self, s, t):
        return self.location.within_boundry(s, t)

    def distance(self, p2):
        return self.location.distance(p2.location)

    def delta(self, p2):
        return self.location.delta(p2.location)

    def initialize_F(self):
        self.ΣF = Vector([0, 0])

    def feel_force(self, F):
        self.ΣF.translate(F)

    def feel_force_from(self, p2, can_merge=False):
        if self not in Particles or p2 not in Particles:
            return False
        # F = GM1M2/R²
        #   pygame.draw.circle(screen, (255, 0, 255), p2.location.get_translate(relative_pos).to_int_list(), int(p2.R))
        R = self.distance(p2)
        if R <= self.R + p2.R and can_merge:
            # create a new particle containing both masses, found at center of mass
            first_index = Particles.index(self)
            sec_index = Particles.index(p2)
            rev_index = max(first_index, sec_index)
            Particles[rev_index] = total_mass([self, p2])
            if rev_index == sec_index:
                Particles.remove(self)
            else:
                Particles.remove(p2)
            return True
        size_of_F = G * self.M * p2.M / (max(R, self.R+p2.R) ** 2)
        direction = p2.delta(self).div(R)  # relative to self. For p2, it's -1 times that. vector with direction of force, size=1
        self.feel_force(direction.multiply(size_of_F))  # Force vector: size of force, direction of force relative to self
        return False


    def move(self):
        a = self.ΣF.div(self.M).get_multiply(time_step)     # calculate acceleration vector
        self.V.translate(a)     # update velocity
        self.location.translate(self.V.get_multiply(time_step))     # move particle
        self.initialize_F()     # clear force(it has already been used)

    def draw(self, surface):
        pygame.draw.circle(surface, self.color, self.location.get_translate(relative_pos).to_int_list(), int(self.R))

    def erase(self, surface):
        pygame.draw.circle(surface, BACKGROUND_COLOR, self.location.get_translate(relative_pos).to_int_list(), int(self.R))

    def momentum(self):
        return self.V.get_multiply(self.M)

    def translate(self, delta):
        self.location.translate(delta)

    def get_translate(self, delta):
        self.location.get_translate(delta)

    def create_copy(self):
        return Particle(self.M, self.P, self.color, self.location.copy(), V=self.V.copy())

def relevant_particles(particles, boundry_1, boundry_2):
    relevant = []
    for p in particles:
        if p.within_boundry(boundry_1, boundry_2):
            relevant.append(p)
    return relevant


def iterate():
    T = total_mass(Particles)
    for f in Forces:
        if f.p1 and f.p2:
            f.enforce()
        else:
            Forces.remove(f)
    for p in Particles:
        p.move()
    return T


def new_iterate():
    global Particles
    T = total_mass(Particles)
    boundry_1 = physics_boundry[0]+T.location
    boundry_2 = physics_boundry[1]+T.location
    Particles = relevant_particles(Particles, boundry_1, boundry_2)
    beg = time.time()
    tree = BarnsHut(Particles, boundry_1, boundry_2)
    #   print(str(tree.nodes()) + "|" + str(tree.depth(0)))
    #   print("time taken to build tree: " + str(time.time()-beg))
    #   tree.depth(0)
    for p in Particles:
        tree.enforce(p)
    for p in Particles:
        p.move()
    return T
    #   return total_mass(Particles)



def take_frame(screen):
    global frame_number
    name = "000000"
    id = str(frame_number)
    name = name[:-len(id)]
    name += id + ".jpg"
    pygame.image.save(screen, Dir + "/" + name)
    frame_number += 1

def total_mass(particles):
    center_of_mass = Vector([0, 0])
    M = 0.0
    momentum = Vector([0, 0])
    color = (0, 0, 0)
    particount = len(particles)
    for p in particles:
        center_of_mass.translate(p.location.get_multiply(p.M))
        M += p.M
        momentum.translate(p.momentum())
        color = [color[i] + p.color[i] for i in range(3)]
    color = [color[i] // particount for i in range(3)]
    center_of_mass.div(M)
    T = Particle(M, P, color, center_of_mass, Moment=momentum)
    return T


def add_particle(M, P, color, location=Vector([0, 0]), V=None, Moment=None):
    new_p = Particle(M, P, color, location, V, Moment)
    add_this_particle(new_p)

def add_this_particle(particle):
    for p in Particles:
        Forces.append(Force(p, particle))
    Particles.append(particle)


def draw_rect(a, b, color):
    g = np.concatenate((a.get_translate(relative_pos).coords.astype(int), b.delta(a).coords.astype(int)))
    pygame.draw.rect(screen, color, g, 1)


def draw_world():
    return
    #   screen.fill((54, 57, 63))
    for p in Particles:
        p.draw(screen)
    #   T.draw(screen)
    #   pygame.display.flip()


def calc_momentum():
    m = Vector([0, 0])
    for p in Particles:
        m.translate(p.momentum())
    return m


def random_place(boundries):
    return Vector([random()*b for b in boundries])


def random_size(boundry):
    return random()*boundry


def range_rand(a, b):
    return (b-a)*random() + a


def random_color():
    return (int(random_size(255)), int(random_size(255)), int(random_size(255)))


def create_random_particle(max_m, P):
    m = random()*max_m
    add_particle(m, P, random_color() , random_place(size), V=Vector([range_rand(-1, 1), range_rand(-1, 1)]))

COLLISION_COEFFICIENT = 0.9
G = 3.0    # gravitational constant
time_step = 1.0
tps = 30000    # ticks per second
frame_rate = 60
wanted_time=10
fps = 60 # frames per second

P = 0.013*16
BACKGROUND_COLOR = (0, 0, 0)#    (54, 57, 63)
relative_pos = Vector([0, 0])
mouse_loc = (0, 0)
Dir = "vid"
gif_name = "result"
frame_number = 0

Forces = []
Particles = []
if not os.path.exists(Dir):
    os.makedirs(Dir)
else:
    for f in os.listdir(Dir):
        os.remove(Dir + "/" + f)

pygame.init()
size = [1200, 1200]
top_left = Vector([0, 0])
size_v = Vector([c for c in size])
screen = pygame.display.set_mode(size)
pygame.display.update()
pygame.display.set_caption("Gravity")


physics_boundry = [size_v*(-10), size_v*11]
# Loop until the user clicks the close button.
done = False
clock = pygame.time.Clock()

screen.fill(BACKGROUND_COLOR)
add_particle(300, P, (255, 255, 255), location=Vector([650, 600]), Moment=Vector([0.5, -1.95]))
for i in range(20):
    create_random_particle(15, P)
#for i in range(10):
#    create_random_particle(20, P)
add_particle(1, P, (255, 255, 0), location=Vector([450, 600]), Moment=Vector([0.0, 1.2]))
add_particle(1, P, (0, 255, 0), location=Vector([250, 600]), Moment=Vector([0.0, 0.75]))
add_particle(1, P, (0, 0, 255), location=Vector([650, 0]), Moment=Vector([-0.5, 0.0]))
print("particles: " + str(len(Particles)) + ", Forces: " + str(len(Forces)))
#   add_particle(10, 1.0, (0, 0, 0), Vector([300, 100]), Moment=Vector([0.2, 0.0]))
beggining = time.time()

prev_momentum = calc_momentum()
while not done:
    #   print("particles: " + str(len(Particles)))
    #   T = new_iterate()
    T = iterate()
    new_momentum = calc_momentum()
    dela_m = new_momentum.delta(prev_momentum)
    """if dela_m.distance(Vector([0, 0])) > 0.001:
        print(dela_m.coords)"""
    prev_momentum = new_momentum
    for p in Particles:
        p.draw(screen)
    pygame.display.flip()
    take_frame(screen)
    """if frame_number ==300:
        done=True"""
    clock.tick(tps)
    #   time.sleep(1.0 / tps)
    #   screen.fill(BACKGROUND_COLOR)
    for p in Particles:
        p.erase(screen)
    #   relative_pos = (size_v.get_multiply(0.5)) + (Particles[0].location*(-1))
    relative_pos = size_v.get_multiply(0.5)+(T.location*(-1))
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            done = True
        else:
            if pygame.mouse.get_pressed()[0]:
                mouse_loc = pygame.mouse.get_pos()
                pos = [((size[i]/2) - mouse_loc[i])/10 for i in range(len(mouse_loc))]
                relative_pos.translate(Vector(pos))
total_time = time.time() - beggining
print(total_time)
img_list = []
imgs_names = sorted(os.listdir(Dir), key=lambda x: int(x[:-4]))

for name in imgs_names:
    img_list.append(Image.open(Dir + "/" + name))

if len(img_list) > 1:
    img_list[0].save(gif_name + ".gif", save_all=True, append_images=img_list[1:], duration=20, loop=0)
elif len(img_list) == 1:
    img_list.save(gif_name + ".gif")
pygame.quit()




# PROBLEM IN CHILD AREA DIVIDING PROBABLY