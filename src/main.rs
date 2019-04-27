#![allow(unused_imports)]
#![feature(euclidean_division)]
#![allow(unused_variables)]

use image;
use nalgebra as na;


use std::f32;
use std::path::Path;
use std::fs::File;

use image::ConvertBuffer;
use na::{Vector3, Vector2, Point3, Unit};
use na::{U1, U2, U3};
use ncollide3d::query::Ray;
use ncollide3d::math::Isometry;
use ncollide3d::shape::Ball;
use ncollide3d::query::ray_internal::ray::RayCast;
use nalgebra::geometry::Isometry3;



const MAX_STEPS: usize = 1024;
const PRECISION: f32 = 0.00001;




#[derive(Copy, Clone, Debug)]
struct Camera {
    origin: Point3<f32>,
    dir: Unit<Vector3<f32>>,
    up: Unit<Vector3<f32>>,
    fov: f32,
    res: (usize, usize),
}

impl Camera {
    fn get_ray(&self, x: usize, y: usize) -> Ray<f32> {
        let x = x as f32;
        let y = y as f32;

        let left = - self.up.cross(&self.dir);
        let dir = *self.dir + *self.up * (y * 2. / self.res.1 as f32 - 1.)
            + left * (x * 2. / self.res.0 as f32 - 1.);

        Ray::new(self.origin, dir.normalize())
    }
}


struct Light {
    position: Point3<f32>,
    intensity: f32,
}
impl Light {
    fn new(position: Point3<f32>, intensity: f32) -> Self {
        Self {
            position,
            intensity,
        }
    }
}

struct Sphere {
    ball: Ball<f32>,
    isometry: Isometry<f32>,
    color: image::Rgb<f32>,
}
impl Sphere {
    fn new(ball: Ball<f32>, isometry: Isometry<f32>, color: image::Rgb<f32>) -> Self {
        Self {
            ball,
            isometry,
            color,
        }
    }
}


struct Scene {
    camera: Camera,
    objects: Vec<Sphere>,
    lights: Vec<Light>,
    background_color: image::Rgb<f32>,
}
impl Scene {
    fn render(&self) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        let mut imgbuf = image::ImageBuffer::new(self.camera.res.0 as u32, self.camera.res.1 as u32);

        let clock = std::time::Instant::now();

        for y in 0..self.camera.res.1 {
            for x in 0..self.camera.res.0 {
                let ray = self.camera.get_ray(x, y);
                let y = self.camera.res.1 - y -1;  // flip vertficaly

                //TODO
                let color = self.trace_ray(&ray);
                *imgbuf.get_pixel_mut(x as u32, y as u32) = color;


                if x < 10 && y < 10 {
                    *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb([250., 250., 250.]);
                }
            }
        }
        let elapsed = clock.elapsed();
        println!("render time: {:?}", elapsed);

        imgbuf
    }

    fn trace_ray(&self, ray: &Ray<f32>) -> image::Rgb<f32> {
        let mut intersections: Vec<_> = self.objects.iter().filter_map(
            |s| s.ball.toi_and_normal_with_ray(&s.isometry, &ray, true).map(|x| (x, s.color))
        ).collect();

        intersections.sort_unstable_by(|(inter1, _color1), (inter2, _color2)| inter1.toi.partial_cmp(&inter2.toi).unwrap());
        match intersections.get(0) {
            Some((inter, color)) => {
                let mut coeff = 0.;
                let intersection_point = ray.origin + ray.dir * inter.toi;

                for light in &self.lights {
                    if self.point_can_see_point(&intersection_point, &light.position) {
                        let vec_to_light = (light.position - intersection_point).normalize();
                        let contribution = &inter.normal.dot(&vec_to_light);
                        if *contribution > 0. {
                            coeff += contribution;
                        }
                    }
                }
                let color = image::Rgb([color[0] * coeff, color[1] * coeff, color[2] * coeff]);
                return color;
            },
            None => return self.background_color,
        }
    }

    fn point_can_see_point(&self, p1: &Point3<f32>, p2: &Point3<f32>) -> bool {
        let ray = Ray::new(p1.clone(), p2 - p1);

        for obj in &self.objects {
            match obj.ball.toi_with_ray(&obj.isometry, &ray, true) {
                None => (),
                Some(t) if t > 0.0001 => return false,
                _ => (),
            }
        }
        true
    }

}  // end impl Scene

fn tone_map(img: &image::ImageBuffer<image::Rgb<f32>, Vec<f32>> ) -> image::ImageBuffer<image::Rgb<u8>, Vec<u8>> {
    let clock = std::time::Instant::now();

    let n = img.width() * img.height();
    
    let mean: f32 = img.pixels().map(|p| p[0].max(p[1]).max(p[2])).sum::<f32>() / n as f32;
    let sqr_mean: f32 = img.pixels().map(|p| p[0].max(p[1]).max(p[2]) * p[0].max(p[1]).max(p[2])).sum::<f32>() / n as f32;

    let variance = sqr_mean - mean * mean;
    let max_exposure = mean + variance.sqrt();


    let mut result: image::ImageBuffer<image::Rgb<u8>, Vec<u8>> = image::ImageBuffer::new(img.width(), img.width());

    for (p_out, p_in) in result.pixels_mut().zip(img.pixels()) {
        for i in 0..3 {
            let ref mut sub_out = p_out[i];
            let sub_in = p_in[i];
            *sub_out = (sub_in.max(256f32)) as u8;
            //*sub_out = (((sub_in as f32 / max_exposure) * 256f32).max(256f32)) as u8;
        }
    }

    let elapsed = clock.elapsed();
    println!("tone mapping time: {:?}", elapsed);

    result
}



fn main() {
    let camera = Camera {
        origin: Point3::origin(),
        dir: Vector3::x_axis(),
        up: Vector3::y_axis(),
        fov: 1.,
        res: (500, 500),
    };

    let objects = vec![
        Sphere::new(Ball::new(1.0), Isometry3::translation(20., 0., 0.), image::Rgb([255., 0., 0.])),
        Sphere::new(Ball::new(2.0), Isometry3::translation(20., 10., 0.), image::Rgb([0., 0., 255.])),
        Sphere::new(Ball::new(8.0), Isometry3::translation(20., 0., 10.), image::Rgb([255., 255., 0.])),
        Sphere::new(Ball::new(20.0), Isometry3::translation(50., 0., 0.), image::Rgb([255., 255., 255.])),
        Sphere::new(Ball::new(10000.0), Isometry3::translation(0., -10020., 0.), image::Rgb([255., 255., 255.])),
    ];

    let lights = vec![
        Light::new(Point3::new(0., 100., 0.), 1.0),
    ];

    let background_color: image::Rgb<f32> = image::Rgb([0., 0., 0.]);

    let scene = Scene {
        camera,
        objects,
        lights,
        background_color,
    };

    let imgbuf = scene.render();
    let imgbuf = tone_map(&imgbuf);
    image::ImageRgb8(imgbuf).save("out.png").unwrap();

}
