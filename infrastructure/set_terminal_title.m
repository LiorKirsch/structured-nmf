    function set_terminal_title(title_string)
    % first check that you are using a terminal 
    % than change the terminal name
    
    if usejava('jvm') && ~feature('ShowFigureWindows')
        system(sprintf('echo -en "\033]0;`hostname` - %s\a"',title_string));
    end

    end