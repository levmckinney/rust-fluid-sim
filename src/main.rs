mod viewer;
use bevy::{prelude::*};

fn main() {
    App::build()
        .add_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_asset::<viewer::MyMaterial>()
        .add_startup_system(viewer::setup.system())
        .add_system(viewer::rotate_camera.system())
        .run();
}