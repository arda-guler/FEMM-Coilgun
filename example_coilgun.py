from autocoilgun import *

def activate_all_forwards(proj, stage):
    z_proj = proj.z_pos
    z_stage = stage.z_pos

    if z_proj > z_stage:
        return 60
    else:
        return 0

def activate_one_forward(proj, stage):
    z_proj = proj.z_pos
    z_stage = stage.z_pos
    maxdist = stage.L * 1.5

    if z_proj - z_stage > 0 and z_stage + maxdist > z_proj:
        return 60
    else:
        return 0

# loader params: idx, r_in, thickness, L, mtl, z_pos
# stage params: idx, N, L, r_in, D_wire, t_ext1, t_ext2, mtl_ext, z_pos
# projectile params: r, L, mtl, z_pos, vel, rho

# this is a 30mm cannon with a muzzle velocity slight above 3m/s
steel_density = 7870
loader = LoadingBlock(0, 16, 75, 100, "1010 Steel", 100)
stage1 = Stage(1, 400, 100, 16, 2.5, 50, 35, "1010 Steel", 0)
stage2 = Stage(2, 400, 100, 16, 2.5, 35, 35, "1010 Steel", -100)
stage3 = Stage(3, 400, 100, 16, 2.5, 35, 25, "1010 Steel", -200)
stage4 = Stage(4, 400, 100, 16, 2.5, 25, 25, "1010 Steel", -300)
stages = [loader, stage1, stage2, stage3, stage4]
gun = Coilgun(stages)
proj = Projectile(15, 100, "1010 Steel", 50, 0, steel_density)
print("Projectile mass:", proj.mass*1000, "g")

PerformAnalysisFull(gun, proj, activate_one_forward, "KST", "KST", 30, 3, 0.0025, 20)
# PerformAnalysisDiscretePosition(gun, proj, activate_one_forward, "KST", "KST", 30, 3, 20)
# PerformAnalysisDiscreteTime(gun, proj, activate_one_forward, "KST", "KST", 30, 3, 0.0025)
