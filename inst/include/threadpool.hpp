/*******************************************************************************
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>

class ThreadPool {
public:
    explicit ThreadPool(size_t thread_count);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
      -> std::future<std::invoke_result_t<F, Args...>>;
    ~ThreadPool();

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop = false;
};

inline ThreadPool::ThreadPool(size_t threads) {
    workers.reserve(threads);
    for (size_t i = 0; i < threads; ++i) {
        workers.emplace_back([this] {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    condition.wait(lock, [this] {
                        return stop || !tasks.empty();
                    });
                    if (stop && tasks.empty())
                        return;
                    task = std::move(tasks.front());
                    tasks.pop();
                }
                task();
            }
        });
    }
}

template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
    -> std::future<std::invoke_result_t<F, Args...>> {
    using return_type = std::invoke_result_t<F, Args...>;

    // wrap the callable+args into a packaged_task
    auto task_ptr = std::make_shared<std::packaged_task<return_type()>>(
        [fn = std::forward<F>(f), ... a = std::forward<Args>(args)]() mutable {
            return fn(a...);
        }
    );

    std::future<return_type> res = task_ptr->get_future();
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        if (stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");
        tasks.emplace([task_ptr]{ (*task_ptr)(); });
    }
    condition.notify_one();
    return res;
}

inline ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (auto &worker : workers)
        worker.join();
}

#endif // THREADPOOL_H_
