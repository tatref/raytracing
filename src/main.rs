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



const MAX_STEPS: usize = 1024;
const PRECISION: f32 = 0.00001;



// TODO: https://ncollide.org/geometric_queries/#ray-casting
#[derive(Copy, Clone, Debug)]
struct Ray {
    org: Vector3<f32>,
    dir: Unit<Vector3<f32>>,
}

impl Ray {
    fn intersect_normal(&self, sphere: &Sphere) -> Option<(f32, Unit<Vector3<f32>>)> {
        let oc = self.org - sphere.center;

        let a = self.dir.dot(&self.dir);
        let b = 2. * oc.dot(&self.dir);
        let c = oc.dot(&oc) - sphere.radius * sphere.radius;

        let delta = b * b - 4. * a * c;

        if delta < 0. {
            None
        }
        else {
            let t = -b -delta.sqrt() / ( 2. * a);
            let hit_point = self.org + t * *self.dir;
            let xx = (hit_point - sphere.center).norm();
            println!("xx: {}", xx);

            let n = Unit::new_normalize(hit_point - sphere.center);
            Some((t, n))
        }
    }
}


#[derive(Copy, Clone, Debug)]
struct Camera {
    org: Vector3<f32>,
    dir: Unit<Vector3<f32>>,
    up: Unit<Vector3<f32>>,
    fov: f32,
    res: (usize, usize),
}

impl Camera {
    fn get_ray(&self, x: usize, y: usize) -> Ray {
        let x = x as f32;
        let y = y as f32;

        let left = - self.up.cross(&self.dir);
        let dir = *self.dir + *self.up * (y * 2. / self.res.1 as f32 - 1.)
            + left * (x * 2. / self.res.0 as f32 - 1.);

        Ray {
            org: self.org,
            dir: Unit::new_normalize(dir),
        }
    }
}


struct Sphere {
    center: Vector3<f32>,
    radius: f32,
    color: image::Rgb<f32>,
}

impl Sphere {
    fn new(center: &Vector3<f32>, radius: f32, color : image::Rgb<f32>) -> Self {
        Self {
            center: center.clone(),
            radius,
            color,
        }
    }
}

struct Plane {
    normal: Unit<Vector3<f32>>,
    dist: f32,
}


struct Scene {
    camera: Camera,
    spheres: Vec<Sphere>,
}

impl Scene {
    fn render(&self) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        let back_color: image::Rgb<f32> = image::Rgb([0., 0., 0.]);

        let mut imgbuf = image::ImageBuffer::new(self.camera.res.0 as u32, self.camera.res.1 as u32);

        let clock = std::time::Instant::now();

        for y in 0..self.camera.res.1 {
            for x in 0..self.camera.res.0 {
                let ray = self.camera.get_ray(x, y);
                let y = self.camera.res.1 - y -1;  // flip vertficaly

                *imgbuf.get_pixel_mut(x as u32, y as u32) = back_color;

                //TODO
                //let color = self.trace_ray(&ray);

                let mut intersections: Vec<_> = self.spheres.iter().filter_map(|s| ray.intersect_normal(&s).map(|(d, n)| (d, n, s.color))).collect();
                intersections.sort_unstable_by(|(dist1, _n1, _color1), (dist2, _n2, _color2)| dist1.partial_cmp(dist2).unwrap());
                match intersections.get(0) {
                    Some((dist, normal, color)) => {
                        let coeff = ray.dir.dot(&normal);
                        let final_color = image::Rgb([color[0] * coeff, color[1] * coeff, color[2] * coeff]);
                        *imgbuf.get_pixel_mut(x as u32, y as u32) = final_color;
                    },
                    None => (),
                }



                if x < 10 && y < 10 {
                    *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb([250., 250., 250.]);
                }
            }
        }
        let elapsed = clock.elapsed();
        println!("render time: {:?}", elapsed);

        imgbuf
    }

}

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
        org: Vector3::zeros(),
        dir: Vector3::x_axis(),
        up: Vector3::y_axis(),
        fov: 1.,
        res: (500, 500),
    };

    let spheres = vec![
        Sphere::new(&Vector3::new(20., 0., 0.), 1.0, image::Rgb([255., 0., 0.])),
        Sphere::new(&Vector3::new(20., 10., 0.), 2.0, image::Rgb([255., 255., 255.])),
        Sphere::new(&Vector3::new(20., 0., 10.), 8.0, image::Rgb([128., 128., 255.])),
    ];
    let scene = Scene {
        camera,
        spheres,
    };

    let imgbuf = scene.render();
    let imgbuf = tone_map(&imgbuf);
    image::ImageRgb8(imgbuf).save("out.png").unwrap();

}
