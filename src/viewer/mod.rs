use bevy::prelude::*;

pub struct RotatingCamera {
    speed: f32,
    distance: f32,
    angle: f32
}

pub struct RotateAround;

/// set up a simple 3D scene
pub fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
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
            mesh: meshes.add(Mesh::from(shape::Cube { size: 1.0 })),
            material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
            transform: Transform::from_translation(Vec3::new(0.0, 1.0, 0.0)),
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
}

pub fn rotate_camera(time: Res<Time>, keyboard_input: Res<Input<KeyCode>>, 
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
    // bound the paddle within the walls
}