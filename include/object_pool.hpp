#ifndef OBJECT_POOL
#define OBJECT_POOL

template <typename T>
class Object_pool {
private:
#ifdef CONFIG_PARALLEL
	tbb::concurrent_queue<T*> pool;
#else
    std::deque<T*> pool;
#endif
public:
    template <typename... Args>
    T* acquire(Args&&... args) {
#ifdef CONFIG_PARALLEL
		T* obj;
		if (!pool.try_pop(obj)) {
			return new T(std::forward<Args>(args)...);
		}
#else
        if (pool.empty()) {
            return new T(std::forward<Args>(args)...);
        }
        T* obj = pool.back();
        pool.pop_back();
#endif
        obj->reset(std::forward<Args>(args)...); // Initialize the object
        return obj;
    }

    void release(T* obj) {
#ifdef CONFIG_PARALLEL
        pool.push(obj);
#else
        pool.push_back(obj);
#endif
    }

    ~Object_pool() {
#ifdef CONFIG_PARALLEL
		T* obj;
		while (pool.try_pop(obj)) {
			delete obj;
		}
#else
        for (T* obj : pool) {
            delete obj;
        }
#endif
    }
};

#endif // !OBJECT_POOL
