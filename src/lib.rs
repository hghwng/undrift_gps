use std::f64::consts::PI;

const PI_X: f64 = std::f64::consts::PI * 3000.0 / 180.0;

fn wgs_encrypt(x: f64, y: f64) -> (f64, f64) {
    let r = 20.0 * (6.0 * PI * y).sin() + 20.0 * (2.0 * PI * y).sin();

    let x_p = 20.0 * (PI * x).sin() + 40.0 * (PI / 3.0 * x).sin();
    let x_q = 160.0 * (PI / 12.0 * x).sin() + 320.0 * (PI / 30.0 * x).sin();
    let x_t = -100.0
        + 2.0 * y
        + 3.0 * x
        + 0.2 * x * x
        + 0.1 * x * y
        + 0.2 * y.abs().sqrt()
        + 2.0 / 3.0 * (r + x_p + x_q);

    let y_p = 20.0 * (PI * y).sin() + 40.0 * (PI / 3.0 * y).sin();
    let y_q = 150.0 * (PI / 12.0 * y).sin() + 300.0 * (PI / 30.0 * y).sin();
    let y_t = 300.0
        + y
        + 2.0 * x
        + 0.1 * y * y
        + 0.1 * x * y
        + 0.1 * y.abs().sqrt()
        + 2.0 / 3.0 * (r + y_p + y_q);

    (x_t, y_t)
}

fn is_outside_china(lat: f64, lon: f64) -> bool {
    lat < 0.8293 || lat > 55.8271 || lon < 72.004 || lon > 137.8347
}

/// Convert a WGS-84 coordinate into GCJ-02
pub fn wgs_to_gcj(lat: f64, lon: f64) -> (f64, f64) {
    /* Krasovsky 1940
       a = 6378245.0, 1/f = 298.3
       b = a * (1 - f)
       ee = (a^2 - b^2) / a^2;
    */
    const A: f64 = 6378245.0;
    const EE: f64 = 0.00669342162296594323;

    if is_outside_china(lat, lon) {
        return (lat, lon);
    }

    let lat_rad = PI / 180.0 * lat;
    let magic = 1.0 - EE * lat_rad.sin().powi(2);

    let (lat_t, lon_t) = wgs_encrypt(lat - 35.0, lon - 105.0);
    let lat_d = (lat_t * 180.0) / (PI * A * (1.0 - EE) / (magic * magic.sqrt()));
    let lon_d = (lon_t * 180.0) / (PI * A / magic.sqrt() * lat_rad.cos());

    (lat + lat_d, lon + lon_d)
}

/// Convert a GCJ-02 coordinate into WGS-84
pub fn gcj_to_wgs(lat: f64, lon: f64) -> (f64, f64) {
    let gcj = (lat, lon);
    let mut wgs = gcj;

    const EPS: f64 = 1e-7;
    const MAX_ROUND: u32 = 10;
    for _ in 0..MAX_ROUND {
        let cur = wgs_to_gcj(wgs.0, wgs.1);
        let delta = ((gcj.0 - cur.0), (gcj.1 - cur.1));
        if delta.0.abs() < EPS && delta.1.abs() < EPS {
            break;
        }

        wgs.0 += delta.0;
        wgs.1 += delta.1;
    }

    wgs
}

/// Convert a GCJ-02 coordinate into BD-09
pub fn gcj_to_bd(lat: f64, lon: f64) -> (f64, f64) {
    let z = (lon * lon + lat * lat).sqrt() + 0.00002 * (PI_X * lat).sin();
    let theta = lat.atan2(lon) + 0.000003 * (PI_X * lon).cos();
    (z * theta.sin() + 0.006, z * theta.cos() + 0.0065)
}

/// Convert a BD-09 coordinate into GCJ-02
pub fn bd_to_gcj(lat: f64, lon: f64) -> (f64, f64) {
    let (lat, lon) = (lat - 0.006, lon - 0.0065);
    let z = (lon * lon + lat * lat).sqrt() - 0.00002 * (PI_X * lat).sin();
    let theta = lat.atan2(lon) - 0.000003 * (PI_X * lon).cos();
    (z * theta.sin(), z * theta.cos())
}

/// Convert a BD-09 coordinate into WGS-84
pub fn bd_to_wgs(lat: f64, lon: f64) -> (f64, f64) {
    let (lat, lon) = bd_to_gcj(lat, lon);
    let (lat, lon) = gcj_to_wgs(lat, lon);
    (lat, lon)
}

/// Convert a WGS-84 coordinate into BD-09
pub fn wgs_to_bd(lat: f64, lon: f64) -> (f64, f64) {
    let (lat, lon) = wgs_to_gcj(lat, lon);
    let (lat, lon) = gcj_to_bd(lat, lon);
    (lat, lon)
}

/// Describes a coordinate system.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeodeticSystem {
    Wgs84,
    Gcj02,
    Bd09,
}

impl GeodeticSystem {
    /// Converts a coordinate to the target system.
    pub fn convert_to(self, target: Self, lat: f64, lon: f64) -> (f64, f64) {
        use GeodeticSystem::*;
        match (self, target) {
            (x, y) if x == y => (lat, lon),
            (Wgs84, Gcj02) => wgs_to_gcj(lat, lon),
            (Wgs84, Bd09) => wgs_to_bd(lat, lon),
            (Gcj02, Wgs84) => gcj_to_wgs(lat, lon),
            (Gcj02, Bd09) => gcj_to_bd(lat, lon),
            (Bd09, Wgs84) => bd_to_wgs(lat, lon),
            (Bd09, Gcj02) => bd_to_gcj(lat, lon),
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    fn loc_assert((x, y): (f64, f64), (p, q): (f64, f64)) {
        let eps_assert = |x: f64, y: f64, var: &str| {
            let diff = (x - y).abs();
            assert!(
                diff < 1e-6,
                "{} = {}, {} expected, diff = {}",
                var,
                x,
                y,
                diff
            );
        };

        eps_assert(x, p, "lat");
        eps_assert(y, q, "lon");
    }

    #[test]
    fn bd_to_gcj() {
        loc_assert(super::bd_to_gcj(39.0, 116.0), (38.994267, 115.993417));
    }

    #[test]
    fn gcj_to_bd() {
        loc_assert(super::gcj_to_bd(39.0, 116.0), (39.005826, 116.006558));
    }

    #[test]
    fn wgs_to_gcj() {
        loc_assert(super::wgs_to_gcj(39.0, 116.0), (39.000886, 116.006018));
    }

    #[test]
    fn gcj_to_wgs() {
        loc_assert(super::gcj_to_wgs(39.0, 116.0), (38.999133, 115.994002));
    }
}
