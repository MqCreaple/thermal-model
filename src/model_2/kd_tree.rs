use std::{ops::Range, slice::from_raw_parts_mut};

use eframe::egui::Vec2;
use ordered_float::NotNan;

use crate::model::{Molecule, MoleculeType};

/// The trait of all objects that has a position (represented by a 2D vector)
/// 
/// The trait is used in KD tree to locate the element
pub trait PositionKey {
    fn get_position(&self) -> Vec2;
}

impl PositionKey for Vec2 {
    fn get_position(&self) -> Vec2 {
        *self
    }
}

impl<T: MoleculeType> PositionKey for Molecule<T> {
    fn get_position(&self) -> Vec2 {
        self.pos
    }
}

/// Inner node of KD tree
/// 
/// The inner node does not contain the x and y bounds of its elements
enum NodeInner<'a, E: PositionKey> {
    Leaf(&'a [E]),
    XNode {
        x_div: f32,
        left: Box<Node<'a, E>>,
        right: Box<Node<'a, E>>,
    },
    YNode {
        y_div: f32,
        up: Box<Node<'a, E>>,
        down: Box<Node<'a, E>>,
    },
}

/// A node of the KD tree
pub struct Node<'a, E: PositionKey> {
    x_bound: Range<f32>,
    y_bound: Range<f32>,
    inner: NodeInner<'a, E>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum KDNodeDivideBy {
    X, Y
}

impl<'a, E: PositionKey + Clone> Node<'a, E> {
    /// Construct a new KD tree from a list of elements
    pub fn new(elems: &'a mut [E], max_size: usize, div_by: KDNodeDivideBy) -> Self {
        // calculate the x and y bounds of all the elements
        let mut x_min = f32::INFINITY;
        let mut x_max = -f32::INFINITY;
        let mut y_min = f32::INFINITY;
        let mut y_max = -f32::INFINITY;
        elems.iter().for_each(|elem| {
            let pos = elem.get_position();
            if pos.x < x_min { x_min = pos.x; }
            if pos.x > x_max { x_max = pos.x; }
            if pos.y < y_min { y_min = pos.y; }
            if pos.y > y_max { y_max = pos.y; }
        });
        let x_bound = x_min..x_max;
        let y_bound = y_min..y_max;
        if elems.len() <= max_size {
            // directly return a leaf node
            Self { x_bound, y_bound, inner: NodeInner::Leaf(elems) }
        } else if div_by == KDNodeDivideBy::X {
            // split the elements by x coordinate
            let len = elems.len();
            // TODO: NaN coordinates are not handled
            let (left, mid, right) = elems.select_nth_unstable_by_key(len / 2, |elem| NotNan::new(elem.get_position().x).unwrap());
            let div = (left[left.len() - 1].get_position().x + mid.get_position().x) / 2.0;
            Self {
                x_bound, y_bound,
                inner: NodeInner::XNode {
                    x_div: div,
                    left: Box::new(Self::new(left, max_size, KDNodeDivideBy::Y)),
                    // TODO: unsafe combination of slices
                    right: Box::new(Self::new(unsafe { from_raw_parts_mut(mid, right.len() + 1) }, max_size, KDNodeDivideBy::Y)),
                },
            }
        } else {
            // split the elements by y coordinate
            let len = elems.len();
            // TODO: NaN coordinates are not handled
            let (up, mid, down) = elems.select_nth_unstable_by_key(len / 2, |elem| NotNan::new(elem.get_position().y).unwrap());
            let div = (up[up.len() - 1].get_position().y + mid.get_position().y) / 2.0;
            Self {
                x_bound, y_bound,
                inner: NodeInner::YNode {
                    y_div: div,
                    up: Box::new(Self::new(up, max_size, KDNodeDivideBy::X)),
                    // TODO: unsafe combination of slices
                    down: Box::new(Self::new(unsafe { from_raw_parts_mut(mid, down.len() + 1) }, max_size, KDNodeDivideBy::X)),
                },
            } 
        }
    }

    /// Find all points in the KD tree that are within a given circle
    pub fn query_circle(&self, center: Vec2, radius: f32) -> Vec<&E> {
        if self.x_bound.start > center.x + radius || self.x_bound.end < center.x - radius ||
            self.y_bound.start > center.y + radius || self.y_bound.end < center.y - radius {
            return Vec::new();
        }
        match &self.inner {
            NodeInner::Leaf(elems) => {
                elems.iter().filter_map(|elem| {
                    if Vec2::length(elem.get_position() - center) <= radius {
                        Some(elem)
                    } else {
                        None
                    }
                }).collect()
            },
            NodeInner::XNode { x_div, left, right } => {
                let overlap_left = center.x - radius < *x_div;
                let overlap_right = center.x + radius > *x_div;
                if overlap_left && overlap_right {
                    let mut result = left.query_circle(center, radius);
                    result.extend(right.query_circle(center, radius));
                    result
                } else if overlap_left {
                    left.query_circle(center, radius)
                } else if overlap_right {
                    right.query_circle(center, radius)
                } else {
                    Vec::new()
                }
            },
            NodeInner::YNode { y_div, up, down } => {
                let overlap_up = center.y - radius < *y_div;
                let overlap_down = center.y + radius > *y_div;
                if overlap_up && overlap_down {
                    let mut result = up.query_circle(center, radius);
                    result.extend(down.query_circle(center, radius));
                    result
                } else if overlap_up {
                    up.query_circle(center, radius)
                } else if overlap_down {
                    down.query_circle(center, radius)
                } else {
                    Vec::new()
                }
            },
        }
    }
}