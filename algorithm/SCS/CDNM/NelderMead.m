% f �Ǻ��������ֻ����һ�� N ά��ʸ����Ϊ��������� ������һ������ֵ
% x0 �� N ά��ʸ���� xerr �� xmin ����Ԫ�صľ������
function [xmin, fmin] = NelderMead(f, x0, xerr)
N = numel(x0); % f �� N Ԫ����
x = zeros(N+1,N); % Ԥ��ֵ
y = zeros(1,N+1);
% ���� N+1 ����ʼ��
x(1,:) = x0;
for ii = 1:N
    x(ii+1,:) = x(1,:);
    if x(1,ii) == 0
        x(ii+1,ii) = 0.00025;
    else
        x(ii+1,ii) = 1.05 * x(1,ii);
    end
end
% ��ѭ��
x_last = x*0;
mask = true(1, N+1); % �ı�Ķ���
while true
    if max(max(abs(x(2:end,:) - x(1,:)))) < xerr % �ж����
        break;
    elseif all(x(:) == x_last(:))
        warning('NelderMead: abs err too small, machine precision reached');
        break;
    else
        x_last = x;
    end
    % ��ֵ������
    for ii = find(mask)
        y(ii) = f(x(ii,:));
    end
    [y, order] = sort(y);
    x = x(order,:);
    m = mean(x(1:N,:)); % ƽ��λ��
    r = 2*m - x(N+1,:); % �����
    f_r = f(r);
    mask(:) = false;
    mask(end) = true;
    if y(1) <= f_r && f_r < y(N) % �� 4 ��
        x(N+1,:) = r; continue;
    elseif f_r < y(1) % �� 5 ��
        s = m + 2*(m - x(N+1,:));
        if f(s) < f_r
            x(N+1,:) = s;
        else
            x(N+1,:) = r;
        end
        continue;
    elseif f_r < y(N+1) % �� 6 ��
        c1 = m + (r - m)*0.5;
        if f(c1) < f_r
            x(N+1,:) = c1; continue;
        end
    else % �� 7 ��
        c2 = m + (x(N+1,:) - m)*0.5;
        if f(c2) < y(N+1)
            x(N+1,:) = c2; continue;
        end
    end
    for jj = 2:N+1 % �� 8 ��
        x(jj,:) = x(1,:) + (x(jj,:) - x(1,:))*0.5;
        mask(jj) = true;
    end
end
% �������
xmin = x(1,:);
fmin = f(xmin);
end
