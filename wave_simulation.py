import numpy as np
from PIL import Image
import time
import pygame



size = (100, 90)


def interaction(left, right, ite):
    try:
        state = pygame.mouse.get_pressed()
        pos = list(pygame.mouse.get_pos())
        scr = 'a'
        if pos[0] > screen_size[0] // 2:
            scr = 'b'
        pos[0] %= screen_size[0] // 2
        pos = (int(pos[0] // shrink_size[0]), int(pos[1] // shrink_size[1]))

        light_a = 1 # 0.5 * (1 + np.cos(ite * np.pi / 10))
        light_b = 1
        dark_a = -1
        dark_b = 0
        if state[0]:
            if scr == 'a':
                left[pos[1]-1: pos[1]+2, pos[0]-1:pos[0]+2] = light_a
            else:
                add_to_arr(left, -right)
                right[pos[1], pos[0]] = light_b
                add_to_arr(left, right)
                """add_to_arr(left, -right)
                right[pos[1] - 1: pos[1] + 2, pos[0] - 1:pos[0] + 2] = light_b
                add_to_arr(left, right)"""

        if state[1]:
            m[pos[1]-1:pos[1]+2, pos[0]-1:pos[0]+2] = wall_fill
            #   m[pos[1]-1: pos[1]+2, pos[0]-1:pos[0]+2] = 1 + (-1 * m[pos[1]-1: pos[1]+2, pos[0]-1:pos[0]+2])
            #   a[pos[1]-1: pos[1]+2, pos[0]-1:pos[0]+2] = 1
            if scr == 'a':
                left[pos[1]-1:pos[1]+2, pos[0]-1:pos[0]+2] = 1-wall_fill
            else:
                add_to_arr(left, -right)
                right[pos[1]-1:pos[1]+2, pos[0]-1:pos[0]+2] = 1
                add_to_arr(left, right)

        if state[2]:
            if scr == 'a':
                left[pos[1]-1: pos[1]+2, pos[0]-1:pos[0]+2] = dark_a
            else:
                add_to_arr(left, -right)
                right[pos[1] - 1: pos[1] + 2, pos[0] - 1:pos[0] + 2] = dark_b
                add_to_arr(left, right)
    except Exception as e:
        print(e)



def initialize(mode='zeros'):
    tmp = (size[0] + 2, size[1] + 2)
    if mode.lower() == 'zeros':
        ar = np.zeros(size)
    if mode.lower() == 'ones':
        ar = np.ones(size)
    elif mode.lower() == 'random':
        ar = np.random.random(size)
    elif mode.lower() == 'peaks':
        ar = np.round(np.random.random(size))
    else:
        ar = np.zeros(size)
    ar[0] = 1
    ar[-1] = 1
    ar[:, 0] = 1
    ar[:, -1] = 1
    return ar


def normalize(a):
    mean = np.mean(a)
    std = max(0.0001, np.std(a))
    return (a-mean)/std


def display(left, right):
    disp_size = (900, 1000)
    dim_a = 1 #  / np.amax(left)
    left = normalize(left)
    img = Image.fromarray((left+1)*127.5*dim_a).resize(disp_size).convert('RGB')
    s = img.size
    mod = img.mode
    img = img.tobytes('raw', 'RGB')
    dim_b = 1.0
    img = pygame.image.frombuffer(img, s, mod)
    screen_a.blit(img, (0, 0))
    img = Image.fromarray(right * 255*dim_b).resize(disp_size).convert('RGB')
    s = img.size
    mod = img.mode
    img = img.tobytes('raw', 'RGB')

    img = pygame.image.frombuffer(img, s, mod)
    screen_a.blit(img, (disp_size[0], 0))
    pygame.display.flip()  # update the display

# interaction

# update function


def update(a, p, m=None):   # Two Dimensional
    """

    :param a: world matrix of heights
    :param p: percentage of height to distribute
    :param m: mask of unavailable areas
    :return:
    """
    delta = np.zeros(a.shape)
    if m is None:   # if no mask is given
        m = np.ones(a.shape)    # mask only the outlines of the map
        m[0] = 0.0
        m[-1] = 0.0
        m[:, 0] = 0.0
        m[:, -1] = 0.0
    for y in range(a.shape[0]):
        for x in range(a.shape[1]):

            pouring = a[y, x] * max(0, min(p, 1)) * m[y, x]   # amount of material to be distributed from this cell
            cells = np.sum(m[y-1: y+2, x-1: x+2])   # amount of cells around(including self) to distribute to

            delta[y, x] -= pouring  # remove the distributing material from cell
            delta[y-1: y+2, x-1: x+2] += m[y-1: y+2, x-1: x+2] * pouring / cells    # distribute material evenly
    a += delta
    return a


def dual_update(a, b, p, m=None):   # Two Dimensional
    """

    :param a: world matrix of heights
    :param p: percentage of height to distribute
    :param m: mask of unavailable areas
    :return:
    """
    delta = np.zeros(a.shape)
    epsil = np.zeros(b.shape)
    if m is None:   # if no mask is given
        m = np.ones(a.shape)    # mask only the outlines of the map
        m[0] = 0.0
        m[-1] = 0.0
        m[:, 0] = 0.0
        m[:, -1] = 0.0
    for y in range(a.shape[0]):
        for x in range(a.shape[1]):


            pouring = a[y, x] * p * m[y, x]   # amount of material to be distributed from this cell
            cells = np.sum(m[y-1: y+2, x-1: x+2])   # amount of cells around(including self) to distribute to

            delta[y, x] -= pouring  # remove the distributing material from cell
            if cells > 0:
                delta[y-1: y+2, x-1: x+2] += m[y-1: y+2, x-1: x+2] * pouring / cells    # distribute material evenly

            pouring = b[y, x] * p * m[y, x]  # amount of material to be distributed from this cell
            cells = np.sum(m[y - 1: y + 2, x - 1: x + 2])  # amount of cells around(including self) to distribute to

            epsil[y, x] -= pouring  # remove the distributing material from cell
            if cells > 0:
                epsil[y - 1: y + 2, x - 1: x + 2] += m[y - 1: y + 2,x - 1: x + 2] * pouring / cells  # distribute material evenly
    a += delta
    b += epsil
    return a, b


def new_wavicle_update(left, right, p, m=None):  # Two Dimensional
    """

    :param a: world matrix of heights
    :param p: percentage of height to distribute
    :param m: mask of unavailable areas
    :return:
    """
    delta = np.zeros(left.shape)
    epsil = np.zeros(right.shape)
    if m is None:  # if no mask is given
        m = np.ones(left.shape)  # mask only the outlines of the map
        m[0] = 0.0
        m[-1] = 0.0
        m[:, 0] = 0.0
        m[:, -1] = 0.0
    for y in range(left.shape[0]):
        for x in range(left.shape[1]):

            pouring = left[y, x] * p * m[y, x]  # amount of material to be distributed from this cell
            cells = np.sum(m[y - 1: y + 2, x - 1: x + 2])  # amount of cells around(including self) to distribute to
            delta[y, x] -= pouring  # remove the distributing material from cell
            if cells > 0:
                delta[y - 1: y + 2, x - 1: x + 2] += m[y - 1: y + 2, x - 1: x + 2] * pouring / cells  # distribute material evenly

    for y in range(right.shape[0]):
        for x in range(right.shape[1]):
            beg = [between(y - 1, 1, left.shape[0] - 4), between(x - 1, 1, left.shape[1] - 4)]
            tmp = left[beg[0]:beg[0] + 3, beg[1]:beg[1] + 3] * (m[beg[0]:beg[0] + 3, beg[1]:beg[1] + 3])
            low = np.unravel_index(np.argmax(tmp), tmp.shape)
            idx = [beg[0] + low[0], beg[1] + low[1]]
            if left[y, x] < left[idx[0], idx[1]]:  # make sure that there actually is a gradient, and not just, if all even, choose first
                epsil[idx[0], idx[1]] += m[y, x] * right[y, x]  # - right[idx[0], idx[1]]
                epsil[y, x] -= m[y, x] * right[y, x]  # - right[idx[0], idx[1]]
                delta[y, x] += a[idx[0], idx[1]] - a[y, x]
                delta[idx[0], idx[1]] += a[y, x] - a[idx[0], idx[1]]
                #   delta[idx[0], idx[1]] += right[y, x]    #    - right[idx[0], idx[1]]
    left += delta + epsil
    right += epsil
    #   left = np.maximum(left, right)
    return left, right


def wavicle_update(left, right, p, m=None):  # Two Dimensional
    """

    :param a: world matrix of heights
    :param p: percentage of height to distribute
    :param m: mask of unavailable areas
    :return:
    """
    delta = np.zeros(left.shape)
    epsil = np.zeros(right.shape)
    if m is None:  # if no mask is given
        m = np.ones(left.shape)  # mask only the outlines of the map
        m[0] = 0.0
        m[-1] = 0.0
        m[:, 0] = 0.0
        m[:, -1] = 0.0
    for y in range(0, left.shape[0]):
        for x in range(0, left.shape[1]):

            pouring = left[y, x] * p * m[y, x]  # amount of material to be distributed from this cell
            cells = np.sum(m[y - 1: y + 2, x - 1: x + 2])  # amount of cells around(including self) to distribute to
            delta[y, x] -= pouring  # remove the distributing material from cell
            if cells > 0:
                delta[y - 1: y + 2, x - 1: x + 2] += m[y - 1: y + 2, x - 1: x + 2] * pouring / cells  # distribute material evenly
            if m[y, x] > 0.0:
                beg = [between(y - 1, 1, left.shape[0] - 4), between(x - 1, 1, left.shape[1] - 4)]
                tmp = left[beg[0]:beg[0] + 3, beg[1]:beg[1] + 3] * (m[beg[0]:beg[0]+3, beg[1]:beg[1]+3])
                low = np.unravel_index(np.argmax(tmp), tmp.shape)
                idx = [beg[0] + low[0], beg[1] + low[1]]
                if left[y, x] < left[idx[0], idx[1]]:  # make sure that there actually is a gradient, and not just, if all even, choose first
                    epsil[idx[0], idx[1]] += right[y, x] #  - right[idx[0], idx[1]]
                    epsil[y, x] -= right[y, x]  # - right[idx[0], idx[1]]
                    #   delta[idx[0], idx[1]] += right[y, x]    #    - right[idx[0], idx[1]]
    left = add_to_arr(left, delta + epsil)
    right = add_to_arr(right, epsil)
    #   left = np.maximum(left, right)
    return left, right


def add_to_arr(ar, delt):
    ar += delt
    return ar



def between(a, low, high):
    return min(max(a, low), high)

a = initialize('zero')
b = initialize('zero')



a += b      # add initial field dents


m = np.ones(a.shape)
m[0] = 0.0
m[-1] = 0.0
m[:, 0] = 0.0
m[:, -1] = 0.0


screen_size = (1800, 1000)
pygame.init()
screen_a = pygame.display.set_mode(screen_size)
c = pygame.time.Clock() # create a clock object for timing

#   fig, ax = plt.subplots(1, 1)
#   im = ax.imshow(np.zeros(a.shape))

shrink_size = (screen_size[0] / (size[1] * 2), screen_size[1] / (size[0]))

fps = 1000.0
interval = 1.0 / fps


ite = 0

wall_fill = 0.0

display_left = a
display_right = b
while True:
    display_left, display_right = new_wavicle_update(display_left, display_right, 1, m)
    display(display_left, display_right)
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            raise SystemExit
    keys = pygame.key.get_pressed()
    if keys[pygame.K_c]:
        wall_fill = 1 - wall_fill
        print('wall_fill: ' + str(wall_fill))
    elif keys[pygame.K_v]:
        wall_fill = abs(wall_fill - 0.5)
        print('wall_fill: ' + str(wall_fill))
    interaction(display_left, display_right, ite)
    time.sleep(interval)
    ite += 1
i = input('finished')



# TO DO: Add second wave Field, where it's values at each spot represent the p of the spot in the original field (Not hard)