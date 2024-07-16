use std::error::Error;
use std::fmt::{Debug, Display};
use std::sync::atomic::{self, AtomicI64};
use std::sync::{Arc, Condvar, Mutex};

pub struct ThreadSyncError(&'static str);

impl Debug for ThreadSyncError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Display for ThreadSyncError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Error for ThreadSyncError {}

/// A simple synchronization tool.
/// 
/// The sync counter is set a given number initially. Then, threads can call
/// `sync()` function to decrease the count by 1 and block the thread until the
/// counter is set to zero.
/// 
/// If any thread panics before calling `sync()`, the counter will be set to -1.
/// In this case, the other threads calling `sync()` will not be blocked, but 
/// their return value will be undefined. Generally, the return value will be an
/// error indicating that another thread panics before calling `sync()`, but it
/// is not guaranteed.
#[derive(Clone)]
pub struct ThreadSync {
    sync: Arc<(AtomicI64, Condvar)>,
    initial_count: i64,
}

impl ThreadSync {
    /// Construct a new thread synchronizer.
    pub fn new(n: u32) -> Self {
        Self {
            sync: Arc::new((AtomicI64::new(n as i64), Condvar::new())),
            initial_count: n as i64,
        }
    }

    /// Blocks the current thread until n threads are waiting.
    /// 
    /// If any thread panics before calling `sync()`, the counter will be set to -1.
    pub fn sync(&self) -> Result<(), ThreadSyncError> {
        // TODO: The code has potential synchronization problem.
        // 1. The counter may be set to -1 before the mutex is locked
        // 2. The last thread might notify all before other threads start to wait for the condition
        let mutex = Mutex::new(());
        let guard = mutex.lock().or(Err(ThreadSyncError("Mutex Guard Poisoned")))?;
        let counter = self.sync.0.fetch_sub(1, atomic::Ordering::AcqRel);
        if counter < 0 {
            return Err(ThreadSyncError("Another thread panics before calling sync()"));
        } else if counter == 1 {
            // reset the counter and notify others
            self.sync.0.store(self.initial_count, atomic::Ordering::Release);
            self.sync.1.notify_all();
        } else {
            // wait for the signal
            let _guard = self.sync.1.wait(guard).or(Err(ThreadSyncError("Error occurred during waiting")))?;
            // _gurad is dropped here
        }
        Ok(())
    }

    

    /// Set a panic hook for the current thread.
    /// 
    /// This is to ensure that a thread that panicked before calling `sync()` will
    /// not block the other threads.
    /// 
    /// `set_panic_hook` should be only called once among all threads, since Rust
    /// does not support setting panic hook for a specific thread. Also, if there
    /// is multiple thread synchronizers, only one of them should call `set_panic_hook`.
    /// 
    /// `set_panic_hook` is a costly function, so avoid it and handle `Result`
    /// manually if possible.
    /// 
    /// # Example
    /// 
    /// ```
    /// use std::thread;
    /// use thermal_model::thread_sync::ThreadSync;
    /// 
    /// let sync = ThreadSync::new(3);
    /// let _guard = sync.set_panic_hook();
    /// let threads = (0..3).map(|i| {
    ///     let sync_clone = sync.clone();
    ///     thread::spawn(move || {
    ///         // Some operations
    ///         if i == 1 {
    ///             panic!("Thread {} panics", i);
    ///         }
    ///         sync_clone.sync().unwrap();
    ///         // Some operations
    ///     })
    /// }).collect::<Vec<_>>();
    /// for thread in threads {
    ///     let result = thread.join();
    ///     println!("{:?}", result);
    /// }
    /// ```
    pub fn set_panic_hook(&self) -> ThreadSyncPanicGuard {
        let sync_clone = self.clone();
        let guard = ThreadSyncPanicGuard {
            old_hook: Some(std::panic::take_hook()),
        };
        std::panic::set_hook(Box::new(move |info| {
            log::error!("Thread panics before synchronization");
            log::error!("detail:");
            log::error!("{}", info);
            sync_clone.sync.0.store(-1, atomic::Ordering::SeqCst);
            sync_clone.sync.1.notify_all();
        }));
        guard
    }
}

pub struct ThreadSyncPanicGuard {
    old_hook: Option<Box<dyn Fn(&std::panic::PanicInfo) + 'static + Send + Sync>>,
}

impl Drop for ThreadSyncPanicGuard {
    fn drop(&mut self) {
        // set the hook back to the old hook
        std::panic::set_hook(std::mem::replace(&mut self.old_hook, None).unwrap());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{thread, time::{Duration, Instant}};

    #[test]
    fn test_sync() {
        let sync = ThreadSync::new(10);
        let threads = (0..10).map(|i| {
            let sync_clone = sync.clone();
            thread::spawn(move || {
                println!("Thread {} initializes", i);
                sync_clone.sync().unwrap();
                println!("Thread {} finishes waiting", i);
            })
        }).collect::<Vec<_>>();
        for th in threads {
            th.join().unwrap();
        }
    }

    #[test]
    fn test_sync_count() {
        const N: u32 = 10;
        let sync = ThreadSync::new(N);
        let clock_before = Instant::now();
        let threads = (0..N).map(|i| {
            let sync_clone = sync.clone();
            thread::spawn(move || {
                println!("Thread {} initializes", i);
                thread::sleep(Duration::from_millis(i as u64 * 20));
                println!("Thread {} finishes sleeping", i);
                sync_clone.sync().unwrap();
                let ret = if (Instant::now() - clock_before).as_millis() < (N - 1) as u128 * 20 {
                    Err(format!("Thread {} is not waiting long enough", i))
                } else {
                    Ok(())
                };
                println!("Thread {} finishes waiting", i);
                ret
            })
        }).collect::<Vec<_>>();
        for th in threads {
            th.join().unwrap().unwrap();
        }
    }

    #[test]
    fn test_sync_panic() {
        // env_logger::builder().filter_level(log::LevelFilter::Debug).init();
        let sync = ThreadSync::new(4);
        let _guard = sync.set_panic_hook();
        let threads = (0..4).map(|i| {
            let sync_clone = sync.clone();
            thread::spawn(move || {
                println!("Thread {} initializes", i);
                thread::sleep(Duration::from_millis(i as u64 * 20));
                if i == 1 {
                    panic!("Thread {} panics", i);
                }
                let _ = sync_clone.sync(); // This may or may not be an `Err`
                println!("Thread {} finishes waiting", i);
            })
        }).collect::<Vec<_>>();
        for (i, th) in threads.into_iter().enumerate() {
            let result = th.join();
            if i == 1 {
                assert!(result.is_err());
            } else {
                assert!(result.is_ok());
            }
        }
    }

    #[test]
    fn test_multiple_sync() {
        let sync = ThreadSync::new(5);
        let threads = (0..5).map(|i| {
            let sync_clone = sync.clone();
            thread::spawn(move || {
                println!("Thread {} initializes", i);
                // first synchronization
                thread::sleep(Duration::from_millis(i as u64 * 20));
                println!("Thread {} finishes sleeping", i);
                sync_clone.sync().unwrap();
                println!("Thread {} syncs", i);
                // second synchronization
                thread::sleep(Duration::from_millis((5 - i as u64) * 20));
                println!("Thread {} finishes sleeping again", i);
                sync_clone.sync().unwrap();
                println!("Thread {} syncs", i);
            })
        }).collect::<Vec<_>>();
        for th in threads {
            th.join().unwrap();
        }
    }
}