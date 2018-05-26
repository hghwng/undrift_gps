# Undrift GPS

[Traversal between the Earth and the Mars](https://en.wikipedia.org/wiki/Restrictions_on_geographic_data_in_China#Coordinate_systems) in just 6 functions! Keep it simple and stupid.

Kudos to [Zili FENG](https://github.com/fengzee-me/)!

```rust
/// Convert a GCJ-02 coordinate into WGS-84
pub fn gcj_to_wgs(lat: f64, lon: f64) -> (f64, f64)

/// Convert a GCJ-02 coordinate into BD-09
pub fn gcj_to_bd(lat: f64, lon: f64) -> (f64, f64)

/// Convert a WGS-84 coordinate into GCJ-02
pub fn wgs_to_gcj(lat: f64, lon: f64) -> (f64, f64)

/// Convert a WGS-84 coordinate into BD-09
pub fn wgs_to_bd(lat: f64, lon: f64) -> (f64, f64)
/// Convert a BD-09 coordinate into GCJ-02
pub fn bd_to_gcj(lat: f64, lon: f64) -> (f64, f64)

/// Convert a BD-09 coordinate into WGS-84
pub fn bd_to_wgs(lat: f64, lon: f64) -> (f64, f64)
```
