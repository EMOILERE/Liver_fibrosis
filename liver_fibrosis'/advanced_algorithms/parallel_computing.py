"""
Parallel Computing Module

Implements advanced parallel computing frameworks for high-performance
biological system modeling including load balancing, distributed memory
management, and communication optimization.

Key Features:
- Multi-threading and multi-processing support
- Recursive coordinate bisection for load balancing
- Distributed memory management with MPI-like interface
- Communication optimization with overlap strategies
- Dynamic task scheduling and work stealing
- Performance monitoring and profiling
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable, Any, Union
from dataclasses import dataclass
from enum import Enum
import threading
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time
import warnings
from queue import Queue, Empty
import psutil


class ParallelizationStrategy(Enum):
    """Types of parallelization strategies."""
    SHARED_MEMORY = "shared_memory"
    DISTRIBUTED_MEMORY = "distributed_memory"
    HYBRID = "hybrid"
    GPU_ACCELERATED = "gpu_accelerated"


class LoadBalancingMethod(Enum):
    """Load balancing methods."""
    STATIC = "static"
    DYNAMIC = "dynamic"
    WORK_STEALING = "work_stealing"
    RECURSIVE_BISECTION = "recursive_bisection"


@dataclass
class ComputationalTask:
    """Container for computational tasks."""
    task_id: str
    function: Callable
    args: Tuple = ()
    kwargs: Dict = None
    priority: int = 0
    estimated_cost: float = 1.0
    dependencies: List[str] = None
    
    def __post_init__(self):
        if self.kwargs is None:
            self.kwargs = {}
        if self.dependencies is None:
            self.dependencies = []


@dataclass
class ProcessorInfo:
    """Information about processor resources."""
    processor_id: int
    core_count: int
    memory_gb: float
    current_load: float = 0.0
    assigned_tasks: List[str] = None
    
    def __post_init__(self):
        if self.assigned_tasks is None:
            self.assigned_tasks = []


class ParallelComputingFramework:
    """
    Main framework for parallel computing in biological simulations.
    """
    
    def __init__(self, strategy: ParallelizationStrategy = ParallelizationStrategy.SHARED_MEMORY,
                 num_workers: Optional[int] = None):
        """
        Initialize parallel computing framework.
        
        Args:
            strategy: Parallelization strategy to use
            num_workers: Number of worker processes/threads
        """
        self.strategy = strategy
        self.num_workers = num_workers or mp.cpu_count()
        
        # Initialize components
        self.load_balancer = LoadBalancer()
        self.memory_manager = DistributedMemoryManager()
        self.comm_optimizer = CommunicationOptimizer()
        self.task_scheduler = TaskScheduler()
        self.performance_monitor = PerformanceMonitor()
        
        # Initialize worker pools
        self.thread_pool = None
        self.process_pool = None
        self._initialize_worker_pools()
        
        # Task management
        self.task_queue = Queue()
        self.result_queue = Queue()
        self.completed_tasks = {}
        
    def _initialize_worker_pools(self):
        """Initialize worker pools based on strategy."""
        if self.strategy in [ParallelizationStrategy.SHARED_MEMORY, ParallelizationStrategy.HYBRID]:
            self.thread_pool = ThreadPoolExecutor(max_workers=self.num_workers)
        
        if self.strategy in [ParallelizationStrategy.DISTRIBUTED_MEMORY, ParallelizationStrategy.HYBRID]:
            self.process_pool = ProcessPoolExecutor(max_workers=self.num_workers)
    
    def submit_task(self, task: ComputationalTask) -> str:
        """
        Submit a computational task for parallel execution.
        
        Args:
            task: Computational task to execute
            
        Returns:
            Task ID for tracking
        """
        # Add task to scheduler
        self.task_scheduler.add_task(task)
        
        # Add to task queue
        self.task_queue.put(task)
        
        return task.task_id
    
    def execute_parallel(self, tasks: List[ComputationalTask],
                        timeout: Optional[float] = None) -> Dict[str, Any]:
        """
        Execute multiple tasks in parallel.
        
        Args:
            tasks: List of computational tasks
            timeout: Maximum execution time
            
        Returns:
            Dictionary of results keyed by task ID
        """
        # Start performance monitoring
        self.performance_monitor.start_monitoring()
        
        # Submit all tasks
        futures = {}
        
        for task in tasks:
            # Check dependencies
            if self._dependencies_satisfied(task):
                future = self._submit_single_task(task)
                futures[task.task_id] = future
            else:
                # Queue for later execution
                self.task_queue.put(task)
        
        # Collect results
        results = {}
        start_time = time.time()
        
        while futures or not self.task_queue.empty():
            # Check for completed futures
            completed_futures = []
            for task_id, future in futures.items():
                if future.done():
                    try:
                        result = future.result()
                        results[task_id] = result
                        self.completed_tasks[task_id] = result
                        completed_futures.append(task_id)
                    except Exception as e:
                        results[task_id] = f"Error: {str(e)}"
                        completed_futures.append(task_id)
            
            # Remove completed futures
            for task_id in completed_futures:
                del futures[task_id]
            
            # Check for tasks with satisfied dependencies
            remaining_tasks = []
            while not self.task_queue.empty():
                try:
                    task = self.task_queue.get_nowait()
                    if self._dependencies_satisfied(task):
                        future = self._submit_single_task(task)
                        futures[task.task_id] = future
                    else:
                        remaining_tasks.append(task)
                except Empty:
                    break
            
            # Re-queue tasks with unsatisfied dependencies
            for task in remaining_tasks:
                self.task_queue.put(task)
            
            # Check timeout
            if timeout and (time.time() - start_time) > timeout:
                warnings.warn("Parallel execution timeout reached")
                break
            
            # Small delay to prevent busy waiting
            time.sleep(0.001)
        
        # Stop performance monitoring
        self.performance_monitor.stop_monitoring()
        
        return results
    
    def _submit_single_task(self, task: ComputationalTask):
        """Submit a single task to appropriate worker pool."""
        if self.strategy == ParallelizationStrategy.SHARED_MEMORY:
            return self.thread_pool.submit(task.function, *task.args, **task.kwargs)
        elif self.strategy == ParallelizationStrategy.DISTRIBUTED_MEMORY:
            return self.process_pool.submit(task.function, *task.args, **task.kwargs)
        elif self.strategy == ParallelizationStrategy.HYBRID:
            # Choose based on task characteristics
            if task.estimated_cost > 10.0:  # CPU-intensive tasks
                return self.process_pool.submit(task.function, *task.args, **task.kwargs)
            else:  # I/O or lightweight tasks
                return self.thread_pool.submit(task.function, *task.args, **task.kwargs)
        else:
            raise ValueError(f"Unsupported parallelization strategy: {self.strategy}")
    
    def _dependencies_satisfied(self, task: ComputationalTask) -> bool:
        """Check if task dependencies are satisfied."""
        for dep_id in task.dependencies:
            if dep_id not in self.completed_tasks:
                return False
        return True
    
    def shutdown(self):
        """Shutdown the parallel computing framework."""
        if self.thread_pool:
            self.thread_pool.shutdown(wait=True)
        if self.process_pool:
            self.process_pool.shutdown(wait=True)


class LoadBalancer:
    """
    Load balancer for distributing computational tasks across processors.
    """
    
    def __init__(self, method: LoadBalancingMethod = LoadBalancingMethod.DYNAMIC):
        """Initialize load balancer."""
        self.method = method
        self.processors = []
        self.task_assignments = {}
        
    def register_processor(self, processor: ProcessorInfo):
        """Register a processor for load balancing."""
        self.processors.append(processor)
    
    def assign_task(self, task: ComputationalTask) -> int:
        """
        Assign a task to the best available processor.
        
        Args:
            task: Task to assign
            
        Returns:
            Processor ID for task assignment
        """
        if self.method == LoadBalancingMethod.STATIC:
            return self._static_assignment(task)
        elif self.method == LoadBalancingMethod.DYNAMIC:
            return self._dynamic_assignment(task)
        elif self.method == LoadBalancingMethod.WORK_STEALING:
            return self._work_stealing_assignment(task)
        elif self.method == LoadBalancingMethod.RECURSIVE_BISECTION:
            return self._recursive_bisection_assignment(task)
        else:
            raise ValueError(f"Unknown load balancing method: {self.method}")
    
    def _static_assignment(self, task: ComputationalTask) -> int:
        """Static load balancing using round-robin."""
        if not self.processors:
            return 0
        
        # Simple round-robin assignment
        processor_id = len(self.task_assignments) % len(self.processors)
        self.task_assignments[task.task_id] = processor_id
        
        return processor_id
    
    def _dynamic_assignment(self, task: ComputationalTask) -> int:
        """Dynamic load balancing based on current processor loads."""
        if not self.processors:
            return 0
        
        # Find processor with minimum load
        min_load = float('inf')
        best_processor = 0
        
        for i, processor in enumerate(self.processors):
            total_load = processor.current_load + sum(
                1.0 for task_id in processor.assigned_tasks 
                if task_id in self.task_assignments
            )
            
            if total_load < min_load:
                min_load = total_load
                best_processor = i
        
        # Assign task
        self.processors[best_processor].assigned_tasks.append(task.task_id)
        self.processors[best_processor].current_load += task.estimated_cost
        self.task_assignments[task.task_id] = best_processor
        
        return best_processor
    
    def _work_stealing_assignment(self, task: ComputationalTask) -> int:
        """Work stealing load balancing."""
        # Initially assign to least loaded processor
        processor_id = self._dynamic_assignment(task)
        
        # Implement work stealing logic (simplified)
        # In practice, this would involve inter-processor communication
        
        return processor_id
    
    def _recursive_bisection_assignment(self, task: ComputationalTask) -> int:
        """Recursive coordinate bisection for spatial decomposition."""
        # This is a simplified implementation
        # In practice, RCB would consider spatial locality of tasks
        
        if not hasattr(task, 'spatial_coordinates'):
            return self._dynamic_assignment(task)
        
        # Use spatial coordinates for assignment
        coords = getattr(task, 'spatial_coordinates', [0, 0])
        
        # Simple bisection based on coordinates
        if len(self.processors) >= 2:
            if coords[0] < 0.5:  # Left half
                processor_id = 0
            else:  # Right half
                processor_id = 1
        else:
            processor_id = 0
        
        self.task_assignments[task.task_id] = processor_id
        return processor_id
    
    def update_processor_load(self, processor_id: int, new_load: float):
        """Update processor load information."""
        if 0 <= processor_id < len(self.processors):
            self.processors[processor_id].current_load = new_load
    
    def get_load_statistics(self) -> Dict[str, float]:
        """Get load balancing statistics."""
        if not self.processors:
            return {}
        
        loads = [p.current_load for p in self.processors]
        
        return {
            'mean_load': np.mean(loads),
            'std_load': np.std(loads),
            'min_load': np.min(loads),
            'max_load': np.max(loads),
            'load_imbalance': (np.max(loads) - np.min(loads)) / (np.mean(loads) + 1e-12)
        }


class DistributedMemoryManager:
    """
    Manager for distributed memory allocation and data sharing.
    """
    
    def __init__(self):
        """Initialize distributed memory manager."""
        self.shared_arrays = {}
        self.memory_pools = {}
        self.allocation_stats = {}
        
    def allocate_shared_array(self, name: str, shape: Tuple[int, ...], 
                            dtype: np.dtype = np.float64) -> np.ndarray:
        """
        Allocate a shared memory array.
        
        Args:
            name: Array identifier
            shape: Array shape
            dtype: Data type
            
        Returns:
            Shared memory array
        """
        try:
            # Use multiprocessing shared memory
            size = np.prod(shape) * np.dtype(dtype).itemsize
            
            # Create shared memory array
            shared_array = mp.Array('d' if dtype == np.float64 else 'f', int(np.prod(shape)))
            
            # Wrap in numpy array
            np_array = np.frombuffer(shared_array.get_obj(), dtype=dtype).reshape(shape)
            
            self.shared_arrays[name] = {
                'array': np_array,
                'shared_obj': shared_array,
                'shape': shape,
                'dtype': dtype
            }
            
            # Update allocation statistics
            self.allocation_stats[name] = {
                'size_bytes': size,
                'allocation_time': time.time()
            }
            
            return np_array
            
        except Exception as e:
            warnings.warn(f"Failed to allocate shared array {name}: {str(e)}")
            # Fallback to regular numpy array
            return np.zeros(shape, dtype=dtype)
    
    def get_shared_array(self, name: str) -> Optional[np.ndarray]:
        """Get a shared memory array by name."""
        if name in self.shared_arrays:
            return self.shared_arrays[name]['array']
        return None
    
    def deallocate_shared_array(self, name: str):
        """Deallocate a shared memory array."""
        if name in self.shared_arrays:
            del self.shared_arrays[name]
            if name in self.allocation_stats:
                del self.allocation_stats[name]
    
    def create_memory_pool(self, pool_name: str, total_size_mb: float):
        """Create a memory pool for efficient allocation."""
        self.memory_pools[pool_name] = {
            'total_size': total_size_mb * 1024 * 1024,  # Convert to bytes
            'allocated': 0,
            'allocations': {}
        }
    
    def get_memory_usage(self) -> Dict[str, float]:
        """Get current memory usage statistics."""
        total_allocated = sum(
            stats['size_bytes'] for stats in self.allocation_stats.values()
        )
        
        # Get system memory info
        memory_info = psutil.virtual_memory()
        
        return {
            'total_allocated_mb': total_allocated / (1024 * 1024),
            'system_total_mb': memory_info.total / (1024 * 1024),
            'system_available_mb': memory_info.available / (1024 * 1024),
            'system_percent_used': memory_info.percent
        }


class CommunicationOptimizer:
    """
    Optimizer for inter-process communication patterns.
    """
    
    def __init__(self):
        """Initialize communication optimizer."""
        self.communication_patterns = {}
        self.message_queues = {}
        self.bandwidth_stats = {}
        
    def setup_communication_pattern(self, pattern_name: str,
                                  sender_ids: List[int],
                                  receiver_ids: List[int]):
        """Setup a communication pattern between processors."""
        self.communication_patterns[pattern_name] = {
            'senders': sender_ids,
            'receivers': receiver_ids,
            'message_count': 0,
            'total_bytes': 0
        }
    
    def send_message(self, pattern_name: str, sender_id: int,
                    receiver_id: int, data: Any):
        """Send a message using optimized communication."""
        if pattern_name not in self.message_queues:
            self.message_queues[pattern_name] = Queue()
        
        message = {
            'sender_id': sender_id,
            'receiver_id': receiver_id,
            'data': data,
            'timestamp': time.time(),
            'size_bytes': self._estimate_message_size(data)
        }
        
        self.message_queues[pattern_name].put(message)
        
        # Update statistics
        if pattern_name in self.communication_patterns:
            self.communication_patterns[pattern_name]['message_count'] += 1
            self.communication_patterns[pattern_name]['total_bytes'] += message['size_bytes']
    
    def receive_message(self, pattern_name: str, receiver_id: int,
                       timeout: Optional[float] = None) -> Optional[Any]:
        """Receive a message from the communication pattern."""
        if pattern_name not in self.message_queues:
            return None
        
        try:
            message = self.message_queues[pattern_name].get(timeout=timeout)
            if message['receiver_id'] == receiver_id:
                return message['data']
            else:
                # Put message back if not for this receiver
                self.message_queues[pattern_name].put(message)
                return None
        except Empty:
            return None
    
    def _estimate_message_size(self, data: Any) -> int:
        """Estimate the size of a message in bytes."""
        if isinstance(data, np.ndarray):
            return data.nbytes
        elif isinstance(data, (list, tuple)):
            return sum(self._estimate_message_size(item) for item in data)
        elif isinstance(data, dict):
            return sum(self._estimate_message_size(v) for v in data.values())
        elif isinstance(data, str):
            return len(data.encode('utf-8'))
        else:
            return 64  # Default estimate
    
    def optimize_communication_schedule(self, pattern_name: str) -> List[Tuple[int, int]]:
        """Optimize communication schedule to minimize contention."""
        if pattern_name not in self.communication_patterns:
            return []
        
        pattern = self.communication_patterns[pattern_name]
        senders = pattern['senders']
        receivers = pattern['receivers']
        
        # Simple optimization: pair senders and receivers to minimize conflicts
        schedule = []
        
        for i, sender in enumerate(senders):
            receiver = receivers[i % len(receivers)]
            schedule.append((sender, receiver))
        
        return schedule


class TaskScheduler:
    """
    Intelligent task scheduler with dependency management.
    """
    
    def __init__(self):
        """Initialize task scheduler."""
        self.tasks = {}
        self.dependency_graph = {}
        self.ready_queue = Queue()
        self.running_tasks = set()
        
    def add_task(self, task: ComputationalTask):
        """Add a task to the scheduler."""
        self.tasks[task.task_id] = task
        
        # Build dependency graph
        self.dependency_graph[task.task_id] = task.dependencies.copy()
        
        # Check if task is ready to run
        if not task.dependencies:
            self.ready_queue.put(task.task_id)
    
    def get_ready_task(self) -> Optional[ComputationalTask]:
        """Get the next ready task for execution."""
        try:
            task_id = self.ready_queue.get_nowait()
            self.running_tasks.add(task_id)
            return self.tasks[task_id]
        except Empty:
            return None
    
    def complete_task(self, task_id: str):
        """Mark a task as completed and update dependencies."""
        if task_id in self.running_tasks:
            self.running_tasks.remove(task_id)
        
        # Check for tasks that can now be scheduled
        for other_task_id, dependencies in self.dependency_graph.items():
            if task_id in dependencies:
                dependencies.remove(task_id)
                
                # If all dependencies satisfied, add to ready queue
                if not dependencies and other_task_id not in self.running_tasks:
                    self.ready_queue.put(other_task_id)
    
    def get_critical_path(self) -> List[str]:
        """Calculate the critical path through the task dependency graph."""
        # Simplified critical path calculation
        # In practice, this would use more sophisticated algorithms
        
        # Find tasks with no dependencies (start nodes)
        start_tasks = [task_id for task_id, deps in self.dependency_graph.items() if not deps]
        
        if not start_tasks:
            return []
        
        # Simple longest path calculation
        longest_path = []
        max_cost = 0
        
        for start_task in start_tasks:
            path, cost = self._calculate_path_cost(start_task)
            if cost > max_cost:
                max_cost = cost
                longest_path = path
        
        return longest_path
    
    def _calculate_path_cost(self, task_id: str, visited: Optional[set] = None) -> Tuple[List[str], float]:
        """Calculate the cost of the longest path from a task."""
        if visited is None:
            visited = set()
        
        if task_id in visited:
            return [], 0.0  # Avoid cycles
        
        visited.add(task_id)
        task = self.tasks[task_id]
        
        # Find dependent tasks
        dependent_tasks = [
            tid for tid, deps in self.dependency_graph.items()
            if task_id in deps
        ]
        
        if not dependent_tasks:
            return [task_id], task.estimated_cost
        
        # Find the longest path among dependents
        max_path = []
        max_cost = 0
        
        for dep_task in dependent_tasks:
            path, cost = self._calculate_path_cost(dep_task, visited.copy())
            if cost > max_cost:
                max_cost = cost
                max_path = path
        
        return [task_id] + max_path, task.estimated_cost + max_cost


class PerformanceMonitor:
    """
    Performance monitor for parallel computing operations.
    """
    
    def __init__(self):
        """Initialize performance monitor."""
        self.monitoring = False
        self.start_time = None
        self.metrics = {}
        self.cpu_usage_history = []
        self.memory_usage_history = []
        
    def start_monitoring(self):
        """Start performance monitoring."""
        self.monitoring = True
        self.start_time = time.time()
        self.metrics = {
            'cpu_percent': [],
            'memory_percent': [],
            'network_io': [],
            'disk_io': []
        }
        
        # Start monitoring thread
        self.monitor_thread = threading.Thread(target=self._monitor_loop)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
    
    def stop_monitoring(self):
        """Stop performance monitoring."""
        self.monitoring = False
        if hasattr(self, 'monitor_thread'):
            self.monitor_thread.join(timeout=1.0)
    
    def _monitor_loop(self):
        """Main monitoring loop."""
        while self.monitoring:
            try:
                # CPU usage
                cpu_percent = psutil.cpu_percent(interval=0.1)
                self.metrics['cpu_percent'].append(cpu_percent)
                
                # Memory usage
                memory_info = psutil.virtual_memory()
                self.metrics['memory_percent'].append(memory_info.percent)
                
                # Network I/O
                net_io = psutil.net_io_counters()
                self.metrics['network_io'].append({
                    'bytes_sent': net_io.bytes_sent,
                    'bytes_recv': net_io.bytes_recv
                })
                
                # Disk I/O
                disk_io = psutil.disk_io_counters()
                if disk_io:
                    self.metrics['disk_io'].append({
                        'read_bytes': disk_io.read_bytes,
                        'write_bytes': disk_io.write_bytes
                    })
                
                time.sleep(0.5)  # Monitor every 0.5 seconds
                
            except Exception as e:
                warnings.warn(f"Performance monitoring error: {str(e)}")
                break
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get performance monitoring summary."""
        if not self.metrics:
            return {}
        
        summary = {}
        
        # CPU statistics
        if self.metrics['cpu_percent']:
            cpu_data = self.metrics['cpu_percent']
            summary['cpu'] = {
                'mean_percent': np.mean(cpu_data),
                'max_percent': np.max(cpu_data),
                'std_percent': np.std(cpu_data)
            }
        
        # Memory statistics
        if self.metrics['memory_percent']:
            memory_data = self.metrics['memory_percent']
            summary['memory'] = {
                'mean_percent': np.mean(memory_data),
                'max_percent': np.max(memory_data),
                'std_percent': np.std(memory_data)
            }
        
        # Execution time
        if self.start_time:
            summary['execution_time'] = time.time() - self.start_time
        
        return summary
    
    def export_metrics(self, filename: str):
        """Export performance metrics to file."""
        import json
        
        export_data = {
            'metrics': self.metrics,
            'summary': self.get_performance_summary(),
            'monitoring_duration': time.time() - self.start_time if self.start_time else 0
        }
        
        with open(filename, 'w') as f:
            json.dump(export_data, f, indent=2, default=str)
