use bevy::{
    prelude::*,
    render::{
        mesh::shape,
        pipeline::{DynamicBinding, PipelineDescriptor, PipelineSpecialization, RenderPipeline},
        render_graph::{base, AssetRenderResourcesNode, RenderGraph},
        renderer::RenderResources,
        shader::{ShaderStage, ShaderStages},
    },
    type_registry::TypeUuid,
};

pub struct RotatingCamera {
    speed: f32,
    distance: f32,
    angle: f32
}

#[derive(RenderResources, Default, TypeUuid)]
#[uuid = "1e08866c-0b8a-437e-8bce-37733b25127e"]
pub struct MyMaterial {
    pub color: Color,
}

const VERTEX_SHADER: &str = r#"
#version 450
layout(location = 0) in vec3 Vertex_Position;
layout(set = 0, binding = 0) uniform Camera {
    mat4 ViewProj;
};
layout(set = 1, binding = 0) uniform Transform {
    mat4 Model;
};
void main() {
    gl_Position = ViewProj * Model * vec4(Vertex_Position, 1.0);
}
"#;

const FRAGMENT_SHADER: &str = r#"
#version 450
layout(location = 0) out vec4 o_Target;

layout(set = 1, binding = 1) uniform MyMaterial_color {
    vec4 color;
};
layout(set = 0, binding = 0) uniform Camera {
    mat4 ViewProj;
};
layout(set = 1, binding = 0) uniform Transform {
    mat4 Model;
};
void main() {
    vec3 posInModel = (inverse(ViewProj*Model)*gl_FragCoord).xyz;
    vec3 forward = normalize(posInModel);
    vec3 s = vec3(-0.1, 0.0, 0.0);
    for(int i = 0; i <= 99; i++) { 
        if (length(s) < 0.05) {
            o_Target = vec4(i*0.01, 1.0, 1.0, 1.0);
            return;
        }
        s = s + forward*0.01;
    }
    o_Target = vec4(1.0, 0.0, 0.0, 1.0);
}
"#;

pub struct RotateAround;

/// set up a simple 3D scene
pub fn setup(
    mut commands: Commands,
    mut pipelines: ResMut<Assets<PipelineDescriptor>>,
    mut shaders: ResMut<Assets<Shader>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<MyMaterial>>,
    mut render_graph: ResMut<RenderGraph>) {

    // Create a new shader pipeline
    let pipeline_handle = pipelines.add(PipelineDescriptor::default_config(ShaderStages {
        vertex: shaders.add(Shader::from_glsl(ShaderStage::Vertex, VERTEX_SHADER)),
        fragment: Some(shaders.add(Shader::from_glsl(ShaderStage::Fragment, FRAGMENT_SHADER))),
    }));

    // Add an AssetRenderResourcesNode to our Render Graph. This will bind MyMaterial resources to our shader
    render_graph.add_system_node(
        "my_material",
        AssetRenderResourcesNode::<MyMaterial>::new(true),
    );

    // Add a Render Graph edge connecting our new "my_material" node to the main pass node. This ensures "my_material" runs before the main pass
    render_graph
        .add_node_edge("my_material", base::node::MAIN_PASS)
        .unwrap();

    // Create a new material
    let material = materials.add(MyMaterial {
        color: Color::rgb(0.0, 0.8, 0.0),
    });

    // add entities to the world
    commands
        // plane
        //.spawn(PbrComponents {
        //    mesh: meshes.add(Mesh::from(shape::Plane { size: 10.0 })),
        //    material: materials.add(Color::rgb(0.3, 0.5, 0.3).into()),
        //    ..Default::default()
        //})
        // cube
        .spawn(MeshComponents {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 1.0 })),
            render_pipelines: RenderPipelines::from_pipelines(vec![RenderPipeline::specialized(
                pipeline_handle,
                // NOTE: in the future you wont need to manually declare dynamic bindings
                PipelineSpecialization {
                    dynamic_bindings: vec![
                        // Transform
                        DynamicBinding {
                            bind_group: 1,
                            binding: 0,
                        },
                        // MyMaterial_color
                        DynamicBinding {
                            bind_group: 1,
                            binding: 1,
                        },
                    ],
                    ..Default::default()
                },
            )]),
            transform: Transform::from_translation(Vec3::new(0.0, 0.0, 0.0)),
            ..Default::default()
        })
        .with(material)
        .with(RotateAround)
        // light
        .spawn(LightComponents {
            transform: Transform::from_translation(Vec3::new(4.0, 8.0, 4.0)),
            ..Default::default()
        })
        // camera
        .spawn(Camera3dComponents {
            transform: Transform::from_translation(Vec3::new(-3.0, 0.0, 0.0))
                .looking_at(Vec3::default(), Vec3::unit_y()),
            ..Default::default()
        }).with(RotatingCamera{speed: 2.0, angle: 0.0, distance: 5.0});
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