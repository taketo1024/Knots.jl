export ReadWriteLock, read_lock, write_lock

mutable struct ReadWriteLock
    readers::Int
    writers::Int
    writing::Bool
    lock::ReentrantLock
    cond::Threads.Condition

    function ReadWriteLock()
        l = ReentrantLock()
        c = Threads.Condition(l)
        new(0, 0, false, l, c)
    end
end

# see https://en.wikipedia.org/wiki/Readers–writer_lock#Using_a_condition_variable_and_a_mutex

function read_lock(f, l::ReadWriteLock)
    lock(l.lock) do
        while l.writers > 0 || l.writing
            wait(l.cond)
        end

        l.readers += 1
    end

    result = f()

    lock(l.lock) do 
        l.readers -= 1
        if l.readers == 0
            notify(l.cond)
        end
    end
    
    result
end

function write_lock(f, l::ReadWriteLock)
    lock(l.lock) do 
        l.writers += 1

        while l.readers > 0 || l.writing
            wait(l.cond)
        end
        
        l.writers -= 1
        l.writing = true
    end

    result = f()

    lock(l.lock) do 
        l.writing = false
        notify(l.cond::Threads.Condition)
    end

    result
end