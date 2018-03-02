function [z,n_k,M,dataset] = read_data(dataset)
datasets = {'text_tomsawyer','text_HOLMES','text_tale','text_WarAndPeace', 'text_MOBY_DICK','text_pride','text_HuckleberryFinn','gene_core','gene_sultan','gene_yang','microdata'};
dataset = datasets{dataset}
switch dataset
    case {'text_tomsawyer','text_HOLMES','text_tale','text_WarAndPeace', 'text_MOBY_DICK','text_pride','text_HuckleberryFinn',}
        
        switch dataset
            case 'text_WarAndPeace'
                %% read "War and Peace, by Leo Tolstoy"
                %'https://www.gutenberg.org/cache/epub/2600/pg2600.txt'
                C = textread('data/pg2600.txt','%s','delimiter','\n');
                for i=1:length(C)
                    if length(C{i})>=3 && strcmp(C{i},'An Anonymous Volunteer, and David Widger')
                        break;
                    end
                end
            case 'text_MOBY_DICK'
                
                %% read "MOBY DICK; OR THE WHALE, By Herman Melville"
                %http://www.gutenberg.org/cache/epub/2701/pg2701.txt
                C = textread('data/pg2701.txt','%s','delimiter','\n');
                for i=1:length(C)
                    if length(C{i})>=3 && strcmp(C{i},'Produced by Daniel Lazarus and Jonesey')
                        break;
                    end
                end
            case 'text_tale'
                C = textread('data/98.txt','%s','delimiter','\n');
                %% read "A TALE OF TWO CITIES, By Charles Dickens"
                for i=1:length(C)
                    if length(C{i})>=3 && strcmp(C{i},'Produced by Judith Boss')
                        break;
                    end
                end
            case 'text_pride'
                C = textread('data/pg1342.txt','%s','delimiter','\n');
                %% read "Pride and Prejudice, by Jane Austen"
                for i=1:length(C)
                    if length(C{i})>=3 && strcmp(C{i},'Produced by Anonymous Volunteers')
                        break;
                    end
                end
            case 'text_tomsawyer'
                C = textread('data/pg74.txt','%s','delimiter','\n');
                %% read "THE ADVENTURES OF TOM SAWYER, By Mark Twain"
                for i=1:length(C)
                    if length(C{i})>=3 && strcmp(C{i},'Produced by David Widger')
                        break;
                    end
                end
            case 'text_HuckleberryFinn'
                C = textread('data/pg76.txt','%s','delimiter','\n');
                %% read "Adventures of Huckleberry Finn, by Mark Twain
                for i=1:length(C)
                    if length(C{i})>=3 && strcmp(C{i},'Produced by David Widger')
                        break;
                    end
                end
            case 'text_HOLMES'
                C = textread('data/pg1661.txt','%s','delimiter','\n');
                %% read "THE ADVENTURES OF SHERLOCK HOLMES, by ARTHUR CONAN DOYLE
                for i=1:length(C)
                    if length(C{i})>=3 && strcmp(C{i},'Produced by an anonymous Project Gutenberg volunteer and Jose Menendez')
                        break;
                    end
                end
        end
        
        for j=i+1:length(C)
            if length(C{j})>=34 && strcmp(C{j}(1:34),'End of the Project Gutenberg EBook')
                break;
            end
        end
        C1=C(i+1:j-1);
        
        
        % C1 = textread('data/dickens.txt','%s','delimiter','\n');
        
        
        C2 = strjoin(C1);
        C3 = regexp(C2, '\w+', 'match');
        C3 = lower(C3);
        C3(cellfun('isempty',C3))=[];
        [Vocabulary,~,z]=unique(C3);
        
        
        
    case {'gene_core','gene_sultan','gene_yang'}
        %http://bowtie-bio.sourceforge.net/recount/countTables/gilad_count_table.txt
        switch dataset
            case 'gene_core'
                fid = fopen('data/core_count_table.txt');
            case 'gene_sultan'
                fid = fopen('data/sultan_count_table.txt');
            case 'gene_yang'
                fid = fopen('data/yang_count_table.txt');
                
        end
        %fid = fopen('data/nagalakshmi_count_table.txt');
        %fid = fopen('data/gilad_count_table.txt');
        
        %fid = fopen('data/yang_count_table.txt');
        %fid = fopen('data/mortazavi_count_table.txt');
        
        
        %fid = fopen('data/katzmouse_count_table.txt');
        %fid = fopen('data/trapnell_count_table.txt');
        A = textscan(fid, '%s','EndOfLine','\r\n','delimiter', ' ');
        A = A{1};
        fclose(fid)
        temp = textscan(A{1},'%s');
        len = length(temp{1});
        X = zeros(length(A)-1,length(temp)-1);
        Gene = cell(length(A)-1,1);
        Label = temp;
        for i=2:length(A)
            temp = textscan(A{i},['%s',repmat('%d ',1,len-1)]);
            for j=2:length(temp)
                X(i-1,j-1) = temp{j}(1);
            end
            Gene(i-1) = temp{1};
        end
        n_k=X(:,1);
        n_k=n_k(n_k>0);
        M = nk_to_m(n_k);
        [z,n_k] = m_i_to_z(M(M>0),find(M));
        [~,~,z]=unique(z);
    case 'microdata'
        m_j = [334 127 63 33 30 20 21 18 16 11 14 13 8 8 4 9 6 9 6 7 9 2 4 3 2 5 6 1 3 3 7 4 1 4 2 1 1 2 2 3 ...
            1 3 2 4 3 1 1 2 1 3 1 2 2 2 1 1 1 1 2 2 3 1 1 2 1 1 1 2 1 1 2 1 1 1 1 1 1 1 2 1 ...
            ones(1,40) 1 1 1 1 1];
        jj =  [1:40 ...
            41:45 47:53 55:58 60 71 73:75 78 79 83 86 88 94 98 107 109 111 112 113 119 131 132 133 144 146 148 ...
            152 158 159 165 166 174 181 183 185 196 206 218 224 244 246 257 272 286 302 310 334 357 378 432 439 505 701 818 836 940 969 1397 1456 1840 2009 2019 2076 2165 3522 3541 ...
            3846 4136 4988 10994 23476];
        M = zeros(1,max(jj));
        M(jj)=m_j;
        [z,n_k] = m_i_to_z(m_j,jj);
        [~,~,z]=unique(z);
end
n_k = full(sparse(1,z,1));
n_k=n_k(n_k>0);
M = nk_to_m(n_k);

