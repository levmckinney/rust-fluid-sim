use bevy::prelude::*;
use nalgebra as na;


fn to_glam(vec: na::Vector3<f32>) -> Vec3 {
    Vec3::new(vec.x, vec.y, vec.z)
}

#[derive(Clone, Copy)]
pub struct SimVolume {
    pub size: na::Vector3<f32>,
    pub center: na::Vector3<f32>,
    pub num_vis_particles: usize
}

pub struct FluidSimViewer {
    pub sim_volume: SimVolume
}

impl Plugin for FluidSimViewer {
    fn build(&self, app: &mut AppBuilder) {
        // add things to your app here
        app.add_resource(self.sim_volume)
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
    sim_volume: ResMut<SimVolume>
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
            mesh: meshes.add(Mesh::from(shape::Cube { size: sim_volume.size.max() })),
            material: materials.add(StandardMaterial {shaded:true, albedo: Color::rgba(1.0,1.0,1.0,0.1) ,..Default::default()}),
            transform: Transform::from_translation(to_glam(sim_volume.center)),
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
            transform: Transform::from_translation(Vec3::new(-3.0, 5.0, 8.0))
                .looking_at(Vec3::default(), Vec3::unit_y()),
            ..Default::default()
        }).with(RotatingCamera{speed: 2.0, angle: 0.0, distance: 10.0});

    for _ in 0..sim_volume.num_vis_particles {
        commands.spawn(
            PbrComponents {
                mesh: meshes.add(Mesh::from(shape::Icosphere { radius: 0.1, subdivisions:1 })),
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

fn update_view_particles(time: Res<Time>, 
                    sim_volume: Res<SimVolume>, 
                    mut particles: Query<(&ViewParticle, &mut Transform)>) {
    
    let width = sim_volume.size.min()/2.0;
    for (_, mut particle_trans) in particles.iter_mut() {
        // TODO update transforms using simulation code
        *particle_trans.translation.x_mut() = (time.seconds_since_startup.cos() as f32)*width;
        *particle_trans.translation.y_mut() = (time.seconds_since_startup.sin() as f32)*width;
        *particle_trans.translation.z_mut() = (time.seconds_since_startup.sin() as f32)*width;
    }
}