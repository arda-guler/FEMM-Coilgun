import femm
import os
import matplotlib.pyplot as plt

PI = 3.14159

class Projectile:
    def __init__(self, r, L, mtl, z_pos=0, vel=0, rho=8000):
        self.r = r
        self.L = L
        self.mtl = mtl
        self.z_pos = z_pos
        self.vel = vel
        self.rho = rho # material density in kg m-3

        self.volume = PI * self.r**2 * (self.L - self.r) # cylindrical segment
        self.volume += 0.5 * (4/3) * PI * self.r**3 # penetrator segment
        self.volume *= 1e-9 # convert from mm3 to m3

        self.mass = self.volume * self.rho # mass in kg

    def get_nodes(self):
        r0 = 0
        r1 = self.r
        z0 = self.z_pos + self.L * 0.5
        z1 = self.z_pos - self.L * 0.5 + self.r
        z2 = self.z_pos - self.L * 0.5

        n0 = [r0, z0]
        n1 = [r1, z0]
        n2 = [r1, z1]
        n3 = [r0, z2]
        nodes = [n0, n1, n2, n3]

        labelpos = [(r0 + r1) * 0.5, (z0 + z1) * 0.5]
        return nodes, labelpos

    def get_accel(self, F):
        return F / self.mass

    def generate_self(self):
        FEMM_GenerateProjectile(self)

class Stage:
    def __init__(self, idx, N, L, r_in, D_wire, t_ext1, t_ext2, mtl_ext, z_pos=0):
        self.idx = idx
        self.N = N
        self.L = L
        self.r_in = r_in
        self.D_wire = D_wire
        self.t_ext1 = t_ext1
        self.t_ext2 = t_ext2
        self.mtl_ext = mtl_ext
        self.z_pos = z_pos

        if D_wire % 1 == 0:
            self.mtl_coil = str(int(D_wire)) + "mm"
        else:
            self.mtl_coil = str(D_wire) + "mm"

        self.r_out = self.r_in + (int(self.N / int(self.L / self.D_wire)) + 1) * self.D_wire # mm, winding outer radius
        self.N_radial = int(self.N / int(self.L / self.D_wire)) + 1 # number of radial layers of winding

        self.I = 0

    def get_nodes(self):
        r0 = self.r_in
        r1 = self.r_out
        r2 = self.r_out + self.t_ext1
        r3 = self.r_out + self.t_ext2

        z0 = self.L * 0.5 + self.z_pos
        z1 = -self.L * 0.5 + self.z_pos

        n0 = [r0, z0]
        n1 = [r1, z0]
        n2 = [r1, z1]
        n3 = [r0, z1]
        n4 = [r2, z0]
        n5 = [r3, z1]
        nodes = [n0, n1, n2, n3, n4, n5]

        coil_label = [(n0[0] + n2[0]) * 0.5, (n0[1] + n2[1]) * 0.5]
        ext_label = [(n1[0] + n5[0]) * 0.5, (n1[1] + n5[1]) * 0.5]

        return nodes, coil_label, ext_label

class LoadingBlock:
    def __init__(self, idx, r_in, thickness, L, mtl, z_pos=0):
        self.idx = idx
        self.r_in = r_in
        self.r_out = r_in + thickness
        self.L = L
        self.mtl = mtl
        self.z_pos = z_pos

    def generate_self(self):
        r0 = self.r_in
        r1 = self.r_out
        z0 = self.z_pos + self.L * 0.5
        z1 = self.z_pos - self.L * 0.5

        n0 = [r0, z0]
        n1 = [r1, z0]
        n2 = [r1, z1]
        n3 = [r0, z1]
        nodes = [n0, n1, n2, n3]

        label_pos = [(r0 + r1) * 0.5, (z0 + z1) * 0.5]

        for n in nodes:
            femm.mi_addnode(n[0], n[1])

        FEMM_ConnectNodes(nodes[0], nodes[1])
        FEMM_ConnectNodes(nodes[1], nodes[2])
        FEMM_ConnectNodes(nodes[2], nodes[3])
        FEMM_ConnectNodes(nodes[3], nodes[0])

        x, y = label_pos[0], label_pos[1]
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop(self.mtl, 1, 1e-3, "<None>", 0, 0, 0)
        femm.mi_clearselected()

class Coilgun:
    def __init__(self, stages):
        self.stages = stages

    def generate_stages(self):
        for idx_s, s in enumerate(self.stages):
            if idx_s == 0:
                s.generate_self()
            else:
                FEMM_GenerateStage(s)

    def set_stage_states(self, statelist):
        for idx_s, s in enumerate(self.stages):
            s.I = statelist[idx_s-1]

def FEMM_Circuit(name, I):
    femm.mi_addcircprop(name, I, 1)

def FEMM_ConnectNodes(n1, n2):
    femm.mi_addsegment(n1[0], n1[1], n2[0], n2[1])

def FEMM_GenerateStage(stage):
    nodes, coil_label_pos, ext_label_pos = stage.get_nodes()

    # geometry
    for n in nodes:
        femm.mi_addnode(n[0], n[1])

    # coil
    FEMM_ConnectNodes(nodes[0], nodes[1])
    FEMM_ConnectNodes(nodes[1], nodes[2])
    FEMM_ConnectNodes(nodes[2], nodes[3])
    FEMM_ConnectNodes(nodes[3], nodes[0])

    # barrel
    FEMM_ConnectNodes(nodes[1], nodes[4])
    FEMM_ConnectNodes(nodes[4], nodes[5])
    FEMM_ConnectNodes(nodes[5], nodes[2])

    # generate circuit
    if stage.I > 0:
        FEMM_Circuit("circuit_" + str(stage.idx), stage.I)

    # coil label
    x, y = coil_label_pos[0], coil_label_pos[1]
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    if stage.I > 0:
        femm.mi_setblockprop(stage.mtl_coil, 1, 1e-3, "circuit_" + str(stage.idx), 0, 0, stage.I)
    else:
        femm.mi_setblockprop(stage.mtl_coil, 1, 1e-3, "<None>", 0, 0, 0)
    femm.mi_clearselected()

    # barrel label
    x, y = ext_label_pos[0], ext_label_pos[1]
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    femm.mi_setblockprop(stage.mtl_ext, 1, 1e-3, "<None>", 0, 0, 0)
    femm.mi_clearselected()

def FEMM_GenerateProjectile(proj):
    nodes, label_pos = proj.get_nodes()

    # geometry
    for n in nodes:
        femm.mi_addnode(n[0], n[1])

    femm.mi_addsegment(nodes[0][0], nodes[0][1], nodes[1][0], nodes[1][1])
    femm.mi_addsegment(nodes[1][0], nodes[1][1], nodes[2][0], nodes[2][1])
    femm.mi_addsegment(nodes[0][0], nodes[0][1], nodes[3][0], nodes[3][1])
    femm.mi_addarc(nodes[3][0], nodes[3][1], nodes[2][0], nodes[2][1], 90, 1)

    # label
    x, y = label_pos[0], label_pos[1]
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    femm.mi_setblockprop(proj.mtl, 1, 1e-3, "<None>", 0, 0, 0)
    femm.mi_clearselected()

    # air label
    x, y = (nodes[0][0] + nodes[1][0]) * 0.5 , nodes[0][1] + 5
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    femm.mi_setblockprop("Air", 1, 1e-3, "<None>", 0, 0, 0)
    femm.mi_clearselected()

def FEMM_GetMaterials(coilgun, projectile):
    mtls = []
    for idx_s, s in enumerate(coilgun.stages):
        if idx_s == 0:
            mtls.append(s.mtl)
        else:
            mtls.append(s.mtl_coil)
            mtls.append(s.mtl_ext)

    mtls.append(projectile.mtl)

    mtls_sanitized = []
    [mtls_sanitized.append(x) for x in mtls if x not in mtls_sanitized]
    mtls_sanitized.append("Air")

    for mtl in mtls_sanitized:
        femm.mi_getmaterial(mtl)

def FEMM_ForceIntegral(proj):
    _, proj_pos = proj.get_nodes()
    femm.mo_selectblock(proj_pos[0], proj_pos[1])
    z_force = femm.mo_blockintegral(19) # in z direction (r should be 0)
    femm.mo_clearblock()
    return z_force

def FEMM_AnalyzeStep(coilgun, state, projectile, dirname, filename, view_padding=30, max_density=1):
    femm.openfemm()
    femm.newdocument(0)
    femm.mi_probdef(0, "millimeters", "axi", 1E-8, 0, 30)

    bitmap_filename = dirname + "/" + filename + ".bmp"
    femm_filename = dirname + "/" + filename + ".fem"
    makedir(dirname)

    FEMM_GetMaterials(coilgun, projectile)

    # create components
    coilgun.set_stage_states(state)
    coilgun.generate_stages()
    projectile.generate_self()

    # create boundary
    femm.mi_makeABC()

    # save and analyze
    femm.mi_saveas(femm_filename)
    femm.mi_createmesh()
    femm.mi_analyse(0)

    # post-process
    femm.mi_loadsolution()
    z_force = FEMM_ForceIntegral(projectile)
    femm.mo_hidegrid()
    femm.mo_showdensityplot(1, 0, max_density, 0, "bmag")

    r_exts = []
    for idx_s, s in enumerate(coilgun.stages):
        if idx_s == 0:
            r_exts.append(s.r_out)
        else:
            r_exts.append(s.r_out + s.t_ext1)
            r_exts.append(s.r_out + s.t_ext2)
    r_max = max(r_exts)

    z_min = coilgun.stages[-1].z_pos - coilgun.stages[-1].L * 0.5
    z_max = coilgun.stages[0].z_pos + coilgun.stages[0].L * 0.5

    x1_zoom = -view_padding
    y1_zoom = z_max + view_padding
    x2_zoom = r_max + view_padding
    y2_zoom = z_min - view_padding
    femm.mo_zoom(x1_zoom, y2_zoom, x2_zoom, y1_zoom)
    femm.mo_savebitmap(bitmap_filename)

    # we are done
    femm.closefemm()

    return z_force, bitmap_filename

def makedir(dirname):
    if not os.path.exists(dirname + "/"):
        os.makedirs(dirname)

def save_movie(files):
    import moviepy.video.io.ImageSequenceClip as mov
    clip = mov.ImageSequenceClip(files, fps=4)
    clip.write_videofile(files[0][:-9] + ".mp4")

def PerformAnalysisDiscreteTime(coilgun, proj, actv_func, dirname, filename, view_padding, max_density, dt=0.01):
    # place projectile at loader
    loader_pos = coilgun.stages[0].z_pos
    proj.z_pos = loader_pos

    muzzle_pos = coilgun.stages[-1].z_pos - coilgun.stages[-1].L * 0.5

    t_forces = []
    t_accels = []
    t_vels = []
    t_bitmaps = []
    times = []
    # discrete time analysis until projectile passes the muzzle
    time = 0
    t_cycle = 0
    dt *= 20
    print("Performing discrete time analysis...")
    while proj.z_pos > muzzle_pos:

        # set stage currents
        currents = []
        for idx_s, s in enumerate(coilgun.stages):
            if not idx_s == 0:
                currents.append(actv_func(proj, s))

        if t_cycle < 10:
            current_filename = filename + "_dt_000" + str(t_cycle)
        elif t_cycle < 100:
            current_filename = filename + "_dt_00" + str(t_cycle)
        elif t_cycle < 1000:
            current_filename = filename + "_dt_0" + str(t_cycle)
        else:
            current_filename = filename + "_dt_" + str(t_cycle)

        force, bitmapname = FEMM_AnalyzeStep(coilgun, currents, proj, dirname, current_filename, view_padding,
                                             max_density)

        accel = force / proj.mass # N / kg = m s-2
        proj.vel += accel * dt # m s-1
        proj.z_pos += proj.vel * dt * 1000 # convert m to mm

        # record data
        t_forces.append(-force)
        t_accels.append(-accel)
        t_vels.append(-proj.vel)
        t_bitmaps.append(bitmapname)
        times.append(time)
        t_cycle += 1
        time += dt
        if t_cycle < 14:
            dt *= 0.794183

    muzzle_vel = -proj.vel

    print("Exporting animation...")
    save_movie(t_bitmaps)

    print("Plotting and exporting results...")
    plt.figure(0)
    plt.plot(times, t_forces)
    plt.xlabel("Time")
    plt.ylabel("Force on Projectile (N)")
    plt.grid()
    plt.savefig(dirname + "/" + filename + "_dt_force.png")

    plt.figure(1)
    plt.plot(times, t_accels)
    plt.xlabel("Time")
    plt.ylabel("Acceleration (m s-2)")
    plt.grid()
    plt.savefig(dirname + "/" + filename + "_dt_accel.png")

    plt.figure(2)
    plt.plot(times, t_vels)
    plt.xlabel("Time")
    plt.ylabel("Velocity (m s-1)")
    plt.grid()
    plt.savefig(dirname + "/" + filename + "_dt_vel.png")

    print("Time spent in barrel:", time, "s")
    print("Muzzle velocity (discrete time analysis):", muzzle_vel, "m/s")
    print("Done.")

def PerformAnalysisDiscretePosition(coilgun, proj, actv_func, dirname, filename, view_padding, max_density, dx=2):
    # place projectile at loader
    loader_pos = coilgun.stages[0].z_pos
    proj.z_pos = loader_pos

    muzzle_pos = coilgun.stages[-1].z_pos - coilgun.stages[-1].L * 0.5

    x_forces = []
    x_accels = []
    x_bitmaps = []
    # discrete position analysis until projectile passes the muzzle
    KE = 0
    x_cycle = 0
    xs = []
    cycles = int((loader_pos - muzzle_pos) / dx) + 1
    print("Performing discrete position analysis...")
    while proj.z_pos > muzzle_pos:

        # set stage currents
        currents = []
        for idx_s, s in enumerate(coilgun.stages):
            if not idx_s == 0:
                currents.append(actv_func(proj, s))

        if x_cycle < 10:
            current_filename = filename + "_dx_000" + str(x_cycle)
        elif x_cycle < 100:
            current_filename = filename + "_dx_00" + str(x_cycle)
        elif x_cycle < 1000:
            current_filename = filename + "_dx_0" + str(x_cycle)
        else:
            current_filename = filename + "_dx_" + str(x_cycle)
        force, bitmapname = FEMM_AnalyzeStep(coilgun, currents, proj, dirname, current_filename, view_padding, max_density)

        accel = force / proj.mass
        proj.z_pos -= dx
        KE += -force * (dx * 1e-3) # convert mm displacement to m

        # record data
        x_forces.append(-force)
        x_accels.append(-accel)
        x_bitmaps.append(bitmapname)
        xs.append(-proj.z_pos - loader_pos)

        x_cycle += 1
        print("Discrete position analysis:", x_cycle / cycles * 100, "%")

    # get muzzle velocity from kinetic energy
    muzzle_vel = (2 * KE / proj.mass)**(0.5)

    print("Exporting animation...")
    save_movie(x_bitmaps)

    print("Plotting and exporting results...")
    plt.figure(0)
    plt.plot(xs, x_forces)
    plt.xlabel("Position Along Barrel (mm)")
    plt.ylabel("Force on Projectile (N)")
    plt.grid()
    plt.savefig(dirname + "/" + filename + "_dx_force.png")

    plt.figure(1)
    plt.plot(xs, x_accels)
    plt.xlabel("Position Along Barrel (mm)")
    plt.ylabel("Projectile Acceleration (m s-2)")
    plt.grid()
    plt.savefig(dirname + "/" + filename + "_dx_accel.png")

    print("Muzzle velocity (discrete position analysis):", muzzle_vel, "m/s")
    print("Done.")

def PerformAnalysisFull(coilgun, proj, actv_func, dirname, filename, view_padding, max_density, dt=0.01, dx=2):
    PerformAnalysisDiscretePosition(coilgun, proj, actv_func, dirname, filename, view_padding, max_density, dx)
    PerformAnalysisDiscreteTime(coilgun, proj, actv_func, dirname, filename, view_padding, max_density, dt)
