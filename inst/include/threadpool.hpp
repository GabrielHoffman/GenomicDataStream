/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/threadpool.hpp
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

// inspired by https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
// replace std::queue< std::function<void()>> with std::queue<std::packaged_task<void()>>
class ThreadPool
{
  public:
    ThreadPool(std::size_t nThreads);
    ~ThreadPool();
    template<class F, class... A>
    decltype(auto) enqueue(F && callable, A &&... arguments); // magic power by modern C++

  private:
    std::vector<std::thread> workers; // keep track of threads
    std::queue<std::packaged_task<void()>> tasks_queue; // use packaged_task instead of function<void()>
    std::mutex mutex_queue;
    std::condition_variable condition; // synchronization
    bool stop;
};

inline ThreadPool::ThreadPool(std::size_t nThreads) : stop(false)
{
    workers.reserve(nThreads);
    for(std::size_t i = 0; i < nThreads; ++i)
    {
        workers.emplace_back(
            [this]
            {
                while(true)
                {
                    std::packaged_task<void()> task; // pack void() func into task
                    {
                        std::unique_lock<std::mutex> lock(mutex_queue);
                        condition.wait(lock, [this] { return stop || !tasks_queue.empty(); });
                        if(tasks_queue.empty() && stop) return;
                        task = std::move(tasks_queue.front()); // front task in the queue moved
                        tasks_queue.pop(); // now pop out the moved front task
                    }
                    task(); // packaged_task created
                }
            });
    }
}

template<class F, class... A>
decltype(auto) ThreadPool::enqueue(F && callable, A &&... arguments)
{
    using ReturnType = std::invoke_result_t<F, A...>;
    std::packaged_task<ReturnType()> task(std::bind(std::forward<F>(callable), std::forward<A>(arguments)...));
    std::future<ReturnType> taskGetFuture = task.get_future();
    {
        std::unique_lock<std::mutex> lock(mutex_queue); // tasks_queue not threading-safe. locking here
        if(stop) // don't allow enqueueing after stopping the pool
            throw std::runtime_error("enqueue on stopped ThreadPool");
        tasks_queue.emplace(std::move(task)); // task moved into queue
    }
    condition.notify_one();
    return taskGetFuture;
}

inline ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(mutex_queue);
        stop = true;
    }
    condition.notify_all();
    for(std::thread & worker : workers) worker.join(); // join all threads
}

#endif // THREADPOOL_H_
