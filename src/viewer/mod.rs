use bevy::prelude::*;
use crate::fluid_sim::{FluidSim, SimConfig};
use crate::fluid_sim::mac_grid::Component;
use crate::fluid_sim::forces::SphereForce;
use nalgebra as na;


fn to_glam(vec: na::Vector3<f64>) -> Vec3 {
    Vec3::new(vec.x as f32, vec.y as f32, vec.z as f32)
}

#[derive(Clone, Copy)]
pub struct SimVolume {
    pub size: na::Vector3<f64>,
    pub center: na::Vector3<f64>,
    pub num_vis_particles: usize
}

pub struct FluidSimViewer {
    pub config: SimConfig
}

impl Plugin for FluidSimViewer {
    fn build(&self, app: &mut AppBuilder) {
        // add things to your app here
        app.add_resource(FluidSim::new(self.config, SphereForce {
                inside_force: na::Vector3::<f64>::new(0.0, 10.0, 0.0),
                outside_force: na::Vector3::<f64>::new(0.0, 0.0, 0.0),
                center: na::Vector3::<f64>::new(0.5, 0.5, 0.5),
                radius: 0.2
            }))
           .add_startup_system(setup.system())
           .add_system(rotate_camera.system())
           .add_system(update_view_particles.system());
    }
}

struct ViewParticle;

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    fluid_sim: Res<FluidSim<SphereForce>>
    ) {
    // add entities to the world
    commands
        // plane
        .spawn(PbrComponents {
            mesh: meshes.add(Mesh::from(shape::Plane { size: 10.0 })),
            material: materials.add(Color::rgb(0.3, 0.5, 0.3).into()),
            ..Default::default()
        })
        // cube
        .spawn(PbrComponents {
            mesh: meshes.add(Mesh::from(shape::Cube { size: (fluid_sim.get_grid().dims.max()/2.0) as f32 })),
            material: materials.add(StandardMaterial {shaded:true, albedo: Color::rgba(1.0,1.0,1.0,0.1) ,..Default::default()}),
            transform: Transform::from_translation(to_glam(fluid_sim.get_grid().dims*0.5)),
            draw: Draw {
                is_transparent: true,
                ..Default::default()
            },
            ..Default::default()
        }).with(RotateAround)
        // light
        .spawn(LightComponents {
            transform: Transform::from_translation(Vec3::new(4.0, 8.0, 4.0)),
            ..Default::default()
        })
        // camera
        .spawn(Camera3dComponents {
            transform: Transform::from_translation(Vec3::new(-3.0, 1.3, 8.0))
                .looking_at(Vec3::default(), Vec3::unit_y()),
            ..Default::default()
        }).with(RotatingCamera{speed: 2.0, angle: 0.0, distance: 1.3});
    // Make view particles 
    for _ in 0..(7*7*7*8) {
        commands.spawn(
            PbrComponents {
                mesh: meshes.add(Mesh::from(shape::Icosphere { radius: 1.0, subdivisions:1 })),
                material: materials.add(StandardMaterial::default()),
                transform: Transform::from_translation(Vec3::zero()),
                ..Default::default()
            }
        ).with(ViewParticle);
    }
}

struct RotatingCamera {
    speed: f32,
    distance: f32,
    angle: f32
}
struct RotateAround;

/// set up a simple 3D scene
fn rotate_camera(time: Res<Time>, 
        keyboard_input: Res<Input<KeyCode>>, 
        mut cam_query: Query<(&mut RotatingCamera, &mut Transform)>, 
        looking_at: Query<(&RotateAround, &Transform)>
    ) {
    let (mut cam, mut cam_transform) = cam_query.iter_mut().last().unwrap();
    let (_, obj_transform) = looking_at.iter().last().unwrap();

    if keyboard_input.pressed(KeyCode::Left) {
        cam.angle += cam.speed*time.delta_seconds;
    }

    if keyboard_input.pressed(KeyCode::Right) {
        cam.angle -= cam.speed*time.delta_seconds;
    }

    let translation = &mut cam_transform.translation;
    *translation.x_mut() = cam.angle.cos()*cam.distance + obj_transform.translation.x();
    *translation.z_mut() = cam.angle.sin()*cam.distance + obj_transform.translation.z();
    *cam_transform = cam_transform.looking_at(obj_transform.translation, Vec3::unit_y());
}

fn update_view_particles(
                    keyboard_input: Res<Input<KeyCode>>,
                    mut fluid_sim: ResMut<FluidSim<SphereForce>>, 
                    mut view_particles: Query<(&ViewParticle, &mut Transform)>) {
    fluid_sim.step();
    let sim_particles = fluid_sim.particles.get_particles();
    for view_particle in view_particles.iter_mut() {
        let (_, mut particle_trans) = view_particle;
        particle_trans.scale = Vec3::zero();
    }
    if keyboard_input.pressed(KeyCode::V) {
        let grad_pressure = fluid_sim.mac_grid.grad_pressure(&fluid_sim.pressure);
        Component::iterator()
            .map(|comp|
                fluid_sim.mac_grid.get_velocity_grid(comp.clone()).cells().map(move |inds| (comp, inds))).flatten()
            .zip(view_particles.iter_mut()).for_each(|((comp, (i,j,k)), view_particle)| {
                let (_, mut particle_trans) = view_particle;
                let position = fluid_sim.mac_grid.u_grid.get_position(i,j,k);
                *particle_trans.translation.x_mut() = position[0] as f32;
                *particle_trans.translation.y_mut() = position[1] as f32;
                *particle_trans.translation.z_mut() = position[2] as f32;
                *particle_trans.scale.x_mut() = 0.01;
                *particle_trans.scale.y_mut() = 0.01;
                *particle_trans.scale.z_mut() = 0.01;
                let vel = (grad_pressure[fluid_sim.mac_grid.get_velocity_ind(comp, i,j,k)] as f32)*0.1;
                match comp {
                    Component::U => {
                        *particle_trans.translation.x_mut() += vel/2.0; 
                        *particle_trans.scale.x_mut() = vel;
                    }
                    Component::V => {
                        *particle_trans.translation.y_mut() += vel/2.0; 
                        *particle_trans.scale.y_mut() = vel;
                    }
                    Component::W => {
                        *particle_trans.translation.z_mut() += vel/2.0; 
                        *particle_trans.scale.z_mut() = vel;
                    }
                }
                 
            });
    } else if keyboard_input.pressed(KeyCode::P) {
        for ((i, j, k), view_particle) in fluid_sim.mac_grid.p_grid.cells().zip(view_particles.iter_mut()) {
            let (_, mut particle_trans) = view_particle;
            let position = fluid_sim.mac_grid.p_grid.get_position(i,j,k);
            *particle_trans.translation.x_mut() = position[0] as f32;
            *particle_trans.translation.y_mut() = position[1] as f32;
            *particle_trans.translation.z_mut() = position[2] as f32;     
            let pressure = fluid_sim.pressure[fluid_sim.mac_grid.get_pressure_ind(i,j,k)]*0.1;
            *particle_trans.scale.x_mut() = pressure as f32;
            *particle_trans.scale.y_mut() = pressure as f32;
            *particle_trans.scale.z_mut() = pressure as f32; 
        }
    } else {
        for (view_particle, sim_particle) in view_particles.iter_mut().zip(sim_particles.into_iter()) {
            let (_, mut particle_trans) = view_particle;
            *particle_trans.translation.x_mut() = sim_particle.position[0] as f32;
            *particle_trans.translation.y_mut() = sim_particle.position[1] as f32;
            *particle_trans.translation.z_mut() = sim_particle.position[2] as f32;  
            *particle_trans.scale.x_mut() = 0.01;
            *particle_trans.scale.y_mut() = 0.01;
            *particle_trans.scale.z_mut() = 0.01;
        }
    }
}